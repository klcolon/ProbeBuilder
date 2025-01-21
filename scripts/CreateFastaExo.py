
import pandas as pd
import subprocess
from pathlib import Path

class CreateBlastDBExo:
    def __init__(self, path, output_path):
        self.path        = path
        self.output_path = Path(output_path)

    def create_fasta_from_exogenous_sequences(self):
            """
            function to create fasta file
            """
            seqs = pd.read_csv(self.path)

            #fill fasta file
            with open(str(self.output_path), "w+") as f:
                for name, seq in seqs.values:
                    f.write(f">{name}\n")
                    f.write(seq + "\n")
            f.close()

            return self

    def create_blast_db(self):
        output = Path(self.output_path.parent) / "exoDB" / "exodb"
        output.parent.mkdir(parents=True, exist_ok=True)
        blast_cmd = f"makeblastdb -in {self.output_path}  -dbtype nucl -parse_seqids -out {output}"
        subprocess.run(blast_cmd, shell=True)

def main(path, output_path):
     CBD = CreateBlastDBExo(path, output_path)
     CBD.create_fasta_from_exogenous_sequences().create_blast_db()

         
#Example
path = "../exogenous_sequences/241101_aavs/blast.csv"
output_path = "../exogenous_sequences/241101_aavs/blast.fa"
main(path, output_path)