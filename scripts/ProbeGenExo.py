#general packages
import pandas     as pd
import numpy      as np
from   Bio.Seq    import Seq
from Bio.SeqUtils import MeltingTemp as mt
from datetime     import datetime
from pathlib      import Path
import os
#for parallel processinng
from   concurrent.futures import ProcessPoolExecutor, as_completed 
#running shell commands
import subprocess

class ProbeGenerator:
    def __init__(self, csv_path, blastdb_rna, blastdb_genome, blast_exogenous, offtarget_overlap=20, probe_size=35, overlap=18, spacing=1):
        self.csv_path        = csv_path        #read in exogenous sequences
        self.probe_size      = probe_size      #length of probe
        self.blastdb_rna     = blastdb_rna     #local blast database for all spliced RNAs
        self.blastdb_genome  = blastdb_genome  #local blast database for genome
        self.blast_exogenous = blast_exogenous #local blast database for exogenous sequences
        self.offtarget_overlap = offtarget_overlap #overlap length for off target
        self.overlap         = overlap         #overlap length for cross hybs within pool
        self.cpu_count       = os.cpu_count() - 1
        self.spacing         = spacing

    def gc_content(self, seq):
        """
        Function to calculate gc content.
        """
        #make sure cases are always upper
        seq = seq.upper()
        #count G's and C's
        G = seq.count("G")
        C = seq.count("C")
        gc = (G+C)/len(seq)
        
        return gc
    
    def generate_fasta_df(self, path = None):
        """
        Read fasta file and convert to Pandas df.
        """
        _list = []
        with open(path, "r") as f:
            temp_gene = []
            temp_seq  = []
            line      = "init"
            while line != "":
                line = f.readline()
                if ">" in line and temp_gene == []:
                    temp_gene.append(line[1:].replace("\n",""))
                elif ">" in line and temp_gene != []:
                    temp_seq = ["".join(temp_seq)]
                    _list.append(temp_gene + temp_seq)
                    temp_gene = []
                    temp_seq  = []
                    temp_gene.append(line[1:].replace("\n",""))
                elif line == "":
                    temp_seq = ["".join(temp_seq)]
                    _list.append(temp_gene + temp_seq)
                else:
                    temp_seq.append(line.replace("\n",""))
            f.close()

        df = pd.DataFrame(_list)
        df.columns = ["Gene","Sequence"]
        
        return df
    
    def sliding_window(self, transcript, kmer=35):
        """
        This function will perform a sliding window of specific kmer size on transcipt sequences
        to generate all possible kmers. These kmers will be filtered based on GC of 0.4-0.60.
        """
        kmer_list = []
        for i in range(len(transcript) - kmer + 1):
            probe = transcript[i:i+kmer]
            probe_gc = self.gc_content(probe)
            if 0.4 <= probe_gc <= 0.65:
                kmer_list.append(probe)
                
        return kmer_list

    def gen_kmers(self, seq):
        window = np.arange(0,(len(seq)-self.overlap+1),1)
        return [seq[i:i+self.overlap] for i in window]

    def kmer_match(self, string1, string2):
        kmer_list = self.gen_kmers(string1)  
        for seq in kmer_list:
            if seq in string2:
                return [[seq, string1, string2]]
        return []

    def check_matches(self, seq_1, seq_2):
        overlap = []
        for i in range(len(seq_1)):
            for j in range(len(seq_2)):
                match = self.kmer_match(seq_1[i], seq_2[j])
                if match:
                    idx = match[0][2].find(match[0][0])
                    overlap.append([match[0][1],match[0][2], idx, self.overlap])
        
        return pd.DataFrame(overlap)

    def check_matches_parallel(self, seq_1, seq_2):
        cores = self.cpu_count
        # Create chunks equivalent to cores
        chunks = np.linspace(0, len(seq_1), cores + 1).astype(int)
        if cores > 1:
            with ProcessPoolExecutor(max_workers=cores) as exe:
                futures = [exe.submit(self.check_matches, seq_1[chunks[i]:chunks[i + 1]], seq_2) for i in range(len(chunks) - 1)]
                results = [fut.result() for fut in as_completed(futures)]
                final_matches = pd.concat(results).reset_index(drop=True)
        else:
            final_matches = self.check_matches(seq_1, seq_2)
        if final_matches.empty:
            print("No cross hybs found...")
        else:
            final_matches.columns = ["Probe", "Off-Target Hit", "Index", "k-mer size"]
            print(f"Cross hybs found: {len(final_matches['Probe'].unique())}")
        return final_matches
        
    def generate_probes_fa(self):
        print("Generating primary probes (ignoring XM and XR)...")
        # Get current time for output directory naming
        current_datetime = datetime.now()
        formatted_datetime = current_datetime.strftime("%Y%m%d%H%M%S")
        home = str(Path.home())
        output_dir = Path(f"{home}/ProbeBuilder/PrimaryProbes/T{formatted_datetime}")
        output_dir.mkdir(parents=True, exist_ok=True)
        self.output_dir = output_dir

        # Read in exogenous sequences
        self.df_seqs = pd.read_csv(self.csv_path)

        # Generate all kmers
        print("Generating all possible probes...")
        gene_kmers_list = []
        for gene, seq in self.df_seqs.values:
            kmers = self.sliding_window(seq, kmer=self.probe_size)
            gene_kmers_list.append([gene, kmers])

        # Write out as fasta file
        with open(f"{output_dir}/primaryblast.fa", "w+") as f:
            for gene in gene_kmers_list:
                for i, seq in enumerate(gene[1]):
                    f.write(f">{gene[0]}_{i}\n")
                    f.write(f"{seq}\n")
        self.blast_fa_location = f"{str(output_dir)}/primaryblast.fa"

        return self

    def blast_primaries_rna(self):
        print("Starting BLAST against RNA database...")
        # BLAST command for RNA database
        blast_cmd = f"blastn -task blastn-short -db {self.blastdb_rna} -query {self.blast_fa_location} -outfmt 6 -out {self.output_dir}/blast_rna_results.txt -word_size 16 -evalue 0.001"
        subprocess.run(blast_cmd, shell=True)
        print("BLAST against RNA database complete...")
        # Read in BLAST output
        try:
            blastout = pd.read_csv(f"{self.output_dir}/blast_rna_results.txt", sep="\t", header=None)
        except:
            blastout = pd.DataFrame()
        
        if len(blastout) == 0:
            self.filtered_fasta_rna = self.generate_fasta_df(self.blast_fa_location)
            # Write the filtered sequences to a new FASTA file for the next BLAST
            filtered_fasta_path = f"{self.output_dir}/filtered_after_rna_blast.fa"
            with open(filtered_fasta_path, "w+") as f:
                for _, row in self.filtered_fasta_rna.iterrows():
                    f.write(f">{row['Gene']}\n")
                    f.write(f"{row['Sequence']}\n")
                f.close()
            # Update the path to the filtered FASTA file
            self.filtered_fasta_location = filtered_fasta_path
            
            return self

        # Filter out off-targets with overlap > 18 nt
        offtarget_idx = []
        for qseqid, sseqid, overlap in blastout[[0, 1, 3]].values:
            if "XM" in sseqid or "XR" in sseqid:
                continue
            else:
                if overlap < self.offtarget_overlap+1:
                    continue
                else:
                    offtarget_idx.append(qseqid)
        print(f"Number of off-target sequences removed after RNA BLAST: {len(set(offtarget_idx))}")

        # Convert FASTA file to pandas dataframe
        fasta = self.generate_fasta_df(self.blast_fa_location)
        # Remove unwanted sequences
        filtered_fasta = fasta[~fasta["Gene"].isin(list(set(offtarget_idx)))].reset_index(drop=True)
        self.final_seq = filtered_fasta

        # Write the filtered sequences to a new FASTA file for the next BLAST
        filtered_fasta_path = f"{self.output_dir}/filtered_after_rna_blast.fa"
        with open(filtered_fasta_path, "w+") as f:
            for _, row in filtered_fasta.iterrows():
                f.write(f">{row['Gene']}\n")
                f.write(f"{row['Sequence']}\n")
            f.close()
        # Update the path to the filtered FASTA file
        self.filtered_fasta_location = filtered_fasta_path

        return self

    def blast_primaries_genome(self):
        print("Starting BLAST against genome database...")
        blast_cmd = f"blastn -task blastn-short -db {self.blastdb_genome} -query {self.filtered_fasta_location} -outfmt 6 -out {self.output_dir}/blast_genome_results.txt -word_size 16 -evalue 0.001"
        subprocess.run(blast_cmd, shell=True)
        print("BLAST against genomic genome database complete...")
        # Read in BLAST output
        try:
             blastout = pd.read_csv(f"{self.output_dir}/blast_genome_results.txt", sep="\t", header=None)
        except:
            blastout = pd.DataFrame()

        if len(blastout) == 0:
            self.filtered_fasta_genomic = self.generate_fasta_df(self.filtered_fasta_location)
            # Write the filtered sequences to a new FASTA file for the next BLAST
            filtered_fasta_path = f"{self.output_dir}/filtered_after_genome_blast.fa"
            with open(filtered_fasta_path, "w+") as f:
                for _, row in self.filtered_fasta_genomic.iterrows():
                    f.write(f">{row['Gene']}\n")
                    f.write(f"{row['Sequence']}\n")
                f.close()
            # Update the path to the filtered FASTA file
            self.filtered_fasta_genomic_location = filtered_fasta_path
            
            return self

        # Filter out off-targets with overlap > 18 nt
        offtarget_idx = []
        for qseqid, sseqid, overlap in blastout[[0, 1, 3]].values:
            # Genomic sequences may have different identifiers, so adjust filtering if necessary
            if sseqid in qseqid:
                continue
            else:
                if overlap < self.offtarget_overlap+1:
                    continue
                else:
                    offtarget_idx.append(qseqid)
        print(f"Number of off-target sequences removed after genomic BLAST: {len(set(offtarget_idx))}")
        # Convert the filtered FASTA file from cDNA BLAST to pandas dataframe
        fasta = self.generate_fasta_df(self.filtered_fasta_location)
        # Remove unwanted sequences based on genomic BLAST results
        self.filtered_fasta_genomic = fasta[~fasta["Gene"].isin(list(set(offtarget_idx)))].reset_index(drop=True)

        # Write the filtered sequences to a new FASTA file for the next BLAST
        filtered_fasta_path = f"{self.output_dir}/filtered_after_genome_blast.fa"
        with open(filtered_fasta_path, "w+") as f:
            for _, row in self.filtered_fasta_genomic.iterrows():
                f.write(f">{row['Gene']}\n")
                f.write(f"{row['Sequence']}\n")
            f.close()
        # Update the path to the filtered FASTA file
        self.filtered_fasta_genomic_location = filtered_fasta_path

        return self

    def blast_primaries_exo(self):
        print("Starting BLAST against exogenous database...")
        # Use the filtered sequences from the genome BLAST as the new query
        blast_cmd = f"blastn -task blastn-short -db {self.blast_exogenous} -query {self.filtered_fasta_genomic_location} -outfmt 6 -out {self.output_dir}/blast_exo_results.txt -word_size 16 -evalue 0.001"
        subprocess.run(blast_cmd, shell=True)
        print("BLAST against exogenous database complete...")
        # Read in BLAST output
        try:
            blastout = pd.read_csv(f"{self.output_dir}/blast_exo_results.txt", sep="\t", header=None)
        except:
            blastout = pd.DataFrame()

        if len(blastout) == 0:
            self.final_seq  = self.generate_fasta_df(self.filtered_fasta_genomic_location)

            return self

        # Filter out off-targets with overlap > 18 nt
        offtarget_idx = []
        for qseqid, sseqid, overlap in blastout[[0, 1, 3]].values:
            # Genomic sequences may have different identifiers, so adjust filtering if necessary
            if sseqid in qseqid:
                continue
            else:
                if overlap < self.offtarget_overlap+1:
                    continue
                else:
                    offtarget_idx.append(qseqid)
        print(f"Number of off-target sequences removed after exogenous BLAST: {len(set(offtarget_idx))}")
        # Convert the filtered FASTA file from genomic BLAST to pandas dataframe
        fasta = self.generate_fasta_df(self.filtered_fasta_genomic_location)
        # Remove unwanted sequences based on genomic BLAST results
        self.final_seq = fasta[~fasta["Gene"].isin(list(set(offtarget_idx)))].reset_index(drop=True)

        return self

    def space_out_probes(self):
        print("Checking probe spacing...")

        # Find start positions of probes in exogenous sequences
        start_list = []
        for name, seq in self.final_seq.values:
            name_prefix = name.split("_")[0]
            exo = self.df_seqs[self.df_seqs.Name == name_prefix].Seqs.iloc[0]
            start = exo.find(seq)
            start_list.append([name, start])

        start_df = pd.DataFrame(start_list, columns=['name', 'start'])
        
        # Remove duplicates and sort
        start_df = start_df.drop_duplicates(subset=['name', 'start'])
        start_df = start_df.sort_values(['start', 'name']).reset_index(drop=True)
        #start_df.to_csv("start.csv")

        final_idx = []
        genes = self.df_seqs["Name"].unique()

        for gene in genes:
            gene_df = start_df[start_df['name'].str.startswith(f"{gene}_")].copy()
            gene_df = gene_df[gene_df['start'] != -1].sort_values('start').reset_index(drop=True)
            if gene_df.empty:
                continue

            # Initialize variables
            prev_end = -1  # End position of the previous probe
            selected_names = []

            for idx, row in gene_df.iterrows():
                name = row['name']
                start = row['start']
                end = start + self.probe_size - 1  # Calculate end position

                if start >= prev_end + self.spacing + 1:
                    # Accept the probe
                    selected_names.append(name)
                    prev_end = end
                else:
                    # Probe overlaps or is too close to the previous one
                    continue

            final_idx.extend(selected_names)

        # Keep the selected probes
        final_gene_seqs = self.final_seq[self.final_seq.Gene.isin(final_idx)].reset_index(drop=True)
        # Remove duplicate sequences
        final_gene_seqs = final_gene_seqs.drop_duplicates(subset="Sequence").reset_index(drop=True)

        # Downsample probes if necessary
        genes = final_gene_seqs.Gene.str.split("_").str[0].unique()
        filtered_df = []
        for gene in genes:
            iso = final_gene_seqs[final_gene_seqs["Gene"].str.contains(f"{gene}_")].reset_index(drop=True)
            if len(iso) > 500:
                filtered_df.append(iso.iloc[:500])
            else:
                filtered_df.append(iso)
            
        filtered_df = pd.concat(filtered_df).reset_index(drop=True)

        self.final_gene_seq = filtered_df

        return self

    def crosshyb_check(self, formamide):
        print("Checking for within-pool cross-hybridizations...")
        # Reverse complement sequences
        rev_comp = [str(Seq(seq).reverse_complement()) for seq in self.final_gene_seq['Sequence']]
        # Check for cross-hybridizations
        matches = self.check_matches_parallel(self.final_gene_seq['Sequence'].values, np.array(rev_comp))
        
        if len(matches) != 0:
            drop_probes = list(set(matches["Probe"].values))
            # Remove probes with cross-hybridizations
            filtered_df = self.final_gene_seq[~self.final_gene_seq['Sequence'].isin(drop_probes)].reset_index(drop=True)
        else:
            filtered_df = self.final_gene_seq.copy()

        # Reverse complement the sequences
        filtered_df["Sequence"] = [str(Seq(seq).reverse_complement()) for seq in filtered_df['Sequence']]
        # Add GC content and melting temperature
        filtered_df["GC"] = filtered_df["Sequence"].apply(self.gc_content)
        # Adjust Tm for formamide concentration
        filtered_df["Tm (formamide adjusted)"] = filtered_df["Sequence"].apply(
            lambda seq: mt.Tm_NN(str(seq), Na=300)) - (0.6 * formamide)

        # Extract the name
        filtered_df["Name"] = filtered_df["Gene"].str.split("_").str[0]

        # Limit the number of probes per gene to a maximum of 50
        max_probes_per_gene = 50
        genes = filtered_df['Name'].unique()
        limited_probes = []
        for gene in genes:
            gene_df = filtered_df[filtered_df['Name'] == gene]
            selected_probes = gene_df.iloc[:max_probes_per_gene]
            limited_probes.append(selected_probes)
        final_df = pd.concat(limited_probes).reset_index(drop=True)

        # Rearrange columns for the final output
        final = final_df[['Name', 'Sequence', 'GC', 'Tm (formamide adjusted)']]

        # Save final probes to CSV
        final.to_csv(str(self.output_dir / "exogenous_probes.csv"), index=False)
        # Save gene info
        gene_info  = final.groupby("Name").size().reset_index()
        gene_info.columns = ["Name", "Number of Probes"]
        gene_info.to_csv(str(self.output_dir / "gene_info.csv"), index=False)

        return self

def main(csv_path, blastdb_rna, blastdb_genome, blast_exogenous, offtarget_overlap =20, probe_size=35, overlap=18, spacing=1, formamide=40):
    # Initialize ProbeGenerator
    PG = ProbeGenerator(csv_path, blastdb_rna, blastdb_genome, blast_exogenous, offtarget_overlap=offtarget_overlap, probe_size=probe_size, overlap=overlap, spacing=spacing)
    PG = PG.generate_probes_fa().blast_primaries_rna().blast_primaries_genome().blast_primaries_exo().space_out_probes()
    PG.crosshyb_check(formamide)
    print("Probes are ready!")

if __name__ == "__main__":
    #example
    csv_path ="../exogenous_sequences/241101_aavs/target.csv"
    blastdb_rna="../mus_musculus/mm39_rna_db/mm39_rna"
    blastdb_genome="../mus_musculus/mm39_genome_db/mm39_genome.fa"
    blast_exogenous="../exogenous_sequences/241101_aavs/exoDB/exodb"
    probe_size=35
    offtarget_overlap=18
    crosshyb_overlap=18
    spacing=0
    formamide=40

    main(csv_path, blastdb_rna, blastdb_genome, blast_exogenous, offtarget_overlap = offtarget_overlap,
        probe_size=probe_size, overlap=crosshyb_overlap, spacing=spacing, formamide=formamide)
    