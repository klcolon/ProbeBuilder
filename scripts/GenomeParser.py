import pandas              as pd
import gffpandas.gffpandas as gffpd
from Bio.Seq               import Seq
from pathlib               import Path
import subprocess

class genome_parser:
    def __init__(self, fasta, gff, output):
        self.fasta    = fasta
        self.gff      = gff
        self.output = output

    def generate_genome_df(self):
        """
        Read fasta file and convert to df.
        """
        _list = []
        with open(self.fasta, "r") as f:
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
            
        genome = pd.DataFrame(_list)
        genome.columns = ["Gene","Sequence"]
        genome = genome[genome["Gene"].str.contains("NC_")].reset_index(drop=True)
        genome["Gene"] = genome["Gene"].str.split(" ").str[0]
        self.genome_df = genome
        
        return self
    
    def generate_gff_df(self):
        """
        Create df from gff file
        """
        gff = gffpd.read_gff3(self.gff).df
        gff = gff[gff.type.isin(["transcript", "gene"])]
        genes = []
        for att in gff.attributes.values:
            start = att.find("gene=")
            sliced = att[start:].split(";")[0].replace("gene=","")
            genes.append(sliced)
        gff["attributes"] = genes
    
        self.gff_df = gff.reset_index(drop=True)
        
        return self
    
    def create_fasta(self):
        """
        function to create fasta file
        """
        genome = self.genome_df
        gff    = self.gff_df
        sequences = []
        for seq_id, start, end, strand, attributes in gff[["seq_id", "start", "end", "strand", "attributes"]].values:
            chrom = genome[genome["Gene"] == seq_id]
            if len(chrom) == 0:
                continue
            gene_seq = chrom["Sequence"].str[start-1:end-1].iloc[0] #-1 for python indexing
            if strand == "-":
                gene_seq = str(Seq(gene_seq).reverse_complement())
                sequences.append([attributes, gene_seq])
            else:
                sequences.append([attributes, gene_seq])
                
        output_path = Path(self.output) / "genome_parsed" /"unspliced_rnas_and_more.fa"
        output_path.parent.mkdir(parents=True, exist_ok=True)
        #fill fasta file
        with open(str(output_path), "w+") as f:
            for i, (name, seq) in enumerate(sequences):
                f.write(f">{name}_{i}\n")
                f.write(seq + "\n")
        f.close()

        self.outputfa = output_path 
    
        return self

    def create_blast_db(self):
        name = Path(self.fasta).parent.name
        output = Path(self.outputfa.parent) / f"{name}_db" / f"{name}_db"
        output.parent.mkdir(parents=True, exist_ok=True)
        blast_cmd = f"makeblastdb -in {self.outputfa}  -dbtype nucl -parse_seqids -out {output}"
        subprocess.run(blast_cmd, shell=True)

def main(input_ref, input_annotations, output_dir):
    gp = genome_parser(input_ref, input_annotations, output_dir)
    gp.generate_genome_df().generate_gff_df().create_fasta().create_blast_db()
        
#example: 
input_ref         = "../mus_musculus/mm39_genomic/GCF_000001635.27_GRCm39_genomic.fna"
input_annotations = "../mus_musculus/mm39_genomic/genomic.gff"
output_dir        = "../mus_musculus"
main(input_ref, input_annotations, output_dir)