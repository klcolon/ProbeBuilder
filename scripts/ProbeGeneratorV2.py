# Author: Katsuya Lex Colon
# Date: 10/25/2024

# Notes: Make sure Gene name matches exactly to human genes including caps. 
# Mouse should be capital in the beginning.
# Note to self: Convert to gui
# Note to self: XM and XR are ignored

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
import multiprocessing    as     mp
#running shell commands
import subprocess

class ProbeGenerator:
    def __init__(self, fasta_path, target, blastdb_cDNA, blastdb_genome, probe_size=35, overlap=18, spacing=1):
        self.fasta_path      = fasta_path   #read in cDNA reference fasta file
        self.probe_size      = probe_size   #length of probe
        self.blastdb_cDNA    = blastdb_cDNA #local blast database for unspliced mRNA
        self.blastdb_genome  = blastdb_genome  #local blast database for genome
        self.overlap         = overlap      #overlap length for cross hybs within pool
        self.target          = target       #RNA targets
        self.cpu_count       = os.cpu_count() - 1
        self.spacing         = spacing

    def generate_fasta_df(self, path = None):
        """
        Read fasta file and convert to Pandas df.
        """
        _list = []
        if path is None:
            path = self.fasta_path
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
    
    def sliding_window(self, transcript, kmer=35):
        """
        This function will perform a sliding window of specific kmer size on transcipt sequences
        to generate all possible kmers. These kmers will be filtered based on GC of 0.4-0.60.
        """
        kmer_list = []
        for i in range(len(transcript) - kmer + 1):
            probe = transcript[i:i+kmer]
            probe_gc = self.gc_content(probe)
            if 0.4 <= probe_gc <= 0.60:
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

        # Convert fasta file to pandas dataframe
        fasta = self.generate_fasta_df()
        # Read in genes of interest
        targets = pd.read_csv(self.target)

        # Find the transcripts and their isoforms, ignoring predicted mRNA isoforms
        sequences = []
        ensemble_ids = []
        genes_arr = targets.values.ravel()
        for gene in genes_arr:
            gene_edit = f"({gene})"
            for _id, seq in fasta.values:
                if gene_edit in _id and "XM_" not in _id and "XR_" not in _id:
                    sequences.append([gene, seq])
                    en_id = _id.split(" ")[0]
                    ensemble_ids.append([gene, en_id])
        self.ensemble_ids = ensemble_ids
        # Convert to dataframe
        df_seqs = pd.DataFrame(sequences)
        self.df_seqs = df_seqs
        # Generate all kmers
        gene_kmers_list = []
        for gene, seq in df_seqs.values:
            kmers = self.sliding_window(seq, kmer=self.probe_size)
            gene_kmers_list.append([gene, kmers])

        # Identify common and variant-specific probes
        print("Looking for common regions for transcript isoforms and unique regions...")
        df_kmers = pd.DataFrame(gene_kmers_list)
        final_kmer_list = []
        for gene in df_kmers[0].unique():
            target_kmers = df_kmers[df_kmers[0] == gene]
            # Check if there are more than one transcript variant
            if len(target_kmers) > 1:
                # Initialize the set of common kmers with the kmers from the first variant
                common_regions = set(target_kmers.iloc[0, 1])
                # Intersect with kmers from other variants
                for i in range(1, len(target_kmers)):
                    common_regions &= set(target_kmers.iloc[i, 1])
                # Add common kmers to the final list
                final_kmer_list.append([f"{gene}_common", list(common_regions)])
                # Add variant-specific kmers for each variant
                for i, seqs in enumerate(target_kmers[1].values):
                    variant_specific_kmers = set(seqs) - common_regions
                    final_kmer_list.append([f"{gene}_{i}", list(variant_specific_kmers)])
            else:
                # Only one variant; all kmers are considered common
                final_kmer_list.append([gene, target_kmers.iloc[0, 1]])

        # Write out as fasta file
        with open(f"{output_dir}/primaryblast.fa", "w+") as f:
            for gene in final_kmer_list:
                for i, seq in enumerate(gene[1]):
                    f.write(f">{gene[0]}_{i}\n")
                    f.write(f"{seq}\n")
        self.blast_fa_location = f"{str(output_dir)}/primaryblast.fa"

        return self

    def blast_primaries_cDNA(self):
        print("Starting BLAST against cDNA database...")
        # BLAST command for cDNA database
        blast_cmd = f"blastn -task blastn-short -db {self.blastdb_cDNA} -query {self.blast_fa_location} -outfmt 6 -out {self.output_dir}/blast_cDNA_results.txt -word_size 16 -evalue 0.001"
        subprocess.run(blast_cmd, shell=True)
        print("BLAST against cDNA database complete...")
        # Read in BLAST output
        blastout = pd.read_csv(f"{self.output_dir}/blast_cDNA_results.txt", sep="\t", header=None)

        # Convert BLAST subject IDs to gene symbols
        gene_sym = []
        for sym_blast in blastout[1].values: 
            found = 0
            for sym, en_id in self.ensemble_ids:
                if sym_blast in en_id:
                    gene_sym.append(sym)
                    found += 1
                    break
            if found == 0:
                gene_sym.append(sym_blast)
        # Reassign names
        blastout[1] = gene_sym
        blastout.to_csv(f"{self.output_dir}/renamed_blast_cDNA_results.txt", sep="\t", index=False, header=False)

        # Filter out off-targets with overlap > 18 nt
        offtarget_idx = []
        for qseqid, sseqid, overlap in blastout[[0, 1, 3]].values:
            if "XM" in sseqid or "XR" in sseqid:
                continue
            if sseqid in qseqid:
                continue
            else:
                if overlap < 19:
                    continue
                else:
                    offtarget_idx.append(qseqid)
        print(f"Number of off-target sequences removed after cDNA BLAST: {len(set(offtarget_idx))}")

        # Convert FASTA file to pandas dataframe
        fasta = self.generate_fasta_df(self.blast_fa_location)
        # Remove unwanted sequences
        filtered_fasta = fasta[~fasta["Gene"].isin(list(set(offtarget_idx)))].reset_index(drop=True)
        self.final_seq_cDNA = filtered_fasta

        # Write the filtered sequences to a new FASTA file for the next BLAST
        filtered_fasta_path = f"{self.output_dir}/filtered_after_cDNA.fa"
        with open(filtered_fasta_path, "w+") as f:
            for idx, row in filtered_fasta.iterrows():
                f.write(f">{row['Gene']}\n")
                f.write(f"{row['Sequence']}\n")
            f.close()
        # Update the path to the filtered FASTA file
        self.filtered_fasta_location = filtered_fasta_path

        return self

    def blast_primaries_genome(self):
        print("Starting BLAST against genomic genome database...")
        # Use the filtered sequences from the cDNA BLAST as the new query
        blast_cmd = f"blastn -task blastn-short -db {self.blastdb_genome} -query {self.filtered_fasta_location} -outfmt 6 -out {self.output_dir}/blast_genomic_genome_results.txt -word_size 16 -evalue 0.001"
        subprocess.run(blast_cmd, shell=True)
        print("BLAST against genomic genome database complete...")
        # Read in BLAST output
        blastout = pd.read_csv(f"{self.output_dir}/blast_genomic_genome_results.txt", sep="\t", header=None)

        # Convert BLAST subject IDs to gene symbols
        gene_sym = []
        for sym_blast in blastout[1].values:
            gene_sym.append(sym_blast.split("_")[0])

        # Reassign names
        blastout[1] = gene_sym
        blastout.to_csv(f"{self.output_dir}/renamed_blast_genome_results.txt", sep="\t", index=False, header=False)

        # Filter out off-targets with overlap > 18 nt
        offtarget_idx = []
        for qseqid, sseqid, overlap in blastout[[0, 1, 3]].values:
            # Genomic sequences may have different identifiers, so adjust filtering if necessary
            if sseqid in qseqid:
                continue
            else:
                if overlap < 19:
                    continue
                else:
                    offtarget_idx.append(qseqid)
        print(f"Number of off-target sequences removed after genomic DNA BLAST: {len(set(offtarget_idx))}")

        # Convert the filtered FASTA file from cDNA BLAST to pandas dataframe
        fasta = self.generate_fasta_df(self.filtered_fasta_location)
        # Remove unwanted sequences based on genomic BLAST results
        self.final_seq = fasta[~fasta["Gene"].isin(list(set(offtarget_idx)))].reset_index(drop=True)

        return self

    def space_out_probes(self):
        print("Checking probe spacing...")

        # Find start positions of probes in transcripts
        start_list = []
        group_df = self.df_seqs.groupby(0)  # Group sequences by gene

        for name, seq in self.final_seq.values:
            name_prefix = name.split("_")[0]

            transcripts = group_df.get_group(name_prefix)[1].values
            # Find the start position in each transcript variant
            for i, transcript in enumerate(transcripts):
                start = transcript.find(seq)
                if start != -1:
                    start_list.append([name, start, f"variant_{i + 1}"])

        start_df = pd.DataFrame(start_list, columns=['name', 'start', 'variant'])
        # Remove duplicates and sort
        start_df = start_df.drop_duplicates(subset=['name', 'start', 'variant'])
        start_df = start_df.sort_values(['start', 'name']).reset_index(drop=True)

        final_idx = []
        genes = self.df_seqs[0].unique()

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
        final_gene_seq = final_gene_seqs.drop_duplicates(subset="Sequence").reset_index(drop=True)

        # Downsample probes if necessary
        genes = final_gene_seq.Gene.str.split("_").str[0].unique()
        filtered_df = []
        for gene in genes:
            iso = final_gene_seq[final_gene_seq["Gene"].str.contains(f"{gene}_")].reset_index(drop=True)
            if len(iso) > 200:
                # Prioritize common probes
                commons = iso[iso["Gene"].str.contains("common")]
                filtered_df.append(commons)
                # Add additional probes up to 200
                no_commons = iso.drop(commons.index).reset_index(drop=True)
                filtered_df.append(no_commons.iloc[:(200 - len(commons))])
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

        # Extract the gene symbol and 'Target Region' information
        filtered_df["Gene_Symbol"] = filtered_df["Gene"].str.split("_").str[0]
        filtered_df["Target Region"] = filtered_df["Gene"].apply(
            lambda x: 'Common' if 'common' in x else 'Variant')

        # Limit the number of probes per gene to a maximum of 50
        max_probes_per_gene = 50
        genes = filtered_df['Gene_Symbol'].unique()
        limited_probes = []
        for gene in genes:
            gene_df = filtered_df[filtered_df['Gene_Symbol'] == gene]
            # Prioritize 'Common' probes
            common_probes = gene_df[gene_df['Target Region'] == 'Common']
            variant_probes = gene_df[gene_df['Target Region'] == 'Variant']
            # Select up to max_probes_per_gene probes
            num_common = min(len(common_probes), max_probes_per_gene)
            selected_probes = common_probes.iloc[:num_common]
            if num_common < max_probes_per_gene:
                num_variant = max_probes_per_gene - num_common
                selected_variant_probes = variant_probes.iloc[:num_variant]
                selected_probes = pd.concat([selected_probes, selected_variant_probes])
            limited_probes.append(selected_probes)
        final_df = pd.concat(limited_probes).reset_index(drop=True)

        # Rearrange columns for the final output
        final = final_df[['Gene_Symbol', 'Sequence', 'Target Region', 'GC', 'Tm (formamide adjusted)']]

        # Save final probes to CSV
        final.to_csv(str(self.output_dir / "probes.csv"), index=False)
        # Save gene info
        final.groupby("Gene_Symbol").size().to_csv(str(self.output_dir / "gene_info.csv"))

        return self

def main(fasta_path, target, blastdb_cDNA, blastdb_genome, probe_size=35, overlap=18, spacing=1, formamide=40):
    # Initialize ProbeGenerator
    PG = ProbeGenerator(fasta_path, target, blastdb_cDNA, blastdb_genome, spacing=spacing, 
                        probe_size=probe_size, overlap=overlap)
    PG = PG.generate_probes_fa().blast_primaries_cDNA().blast_primaries_genome().space_out_probes()
    PG.crosshyb_check(formamide)
    print("Probes are ready!")


if __name__ == "__main__":
    #example
    path_to_full_length_transcriptome = "homo_sapiens/human_db/GRCh38.p14.rna.fna"
    csv_of_mRNA_targets               = "genes_files/2genes.csv"
    path_to_blastdb_RNA               = "homo_sapiens/hg38_cDNA/hg38_cDNA"
    path_to_blastdb_Genomic           = "homo_sapiens/hg38_genome/grch38_genome"
    probe_length                      = 35
    overlap_between_probes            = 18
    spacing_between_probes            = 1
    formamide_use                     = 45

    main(path_to_full_length_transcriptome, csv_of_mRNA_targets, path_to_blastdb_RNA,
        path_to_blastdb_Genomic , probe_length, overlap_between_probes, 
        spacing_between_probes , formamide_use)