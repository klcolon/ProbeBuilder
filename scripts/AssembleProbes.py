import pandas as pd
import numpy as np
from pathlib import Path
##edit

class AssembleProbes:
    def __init__(self, primary, codebook, readouts):
        self.primary = pd.read_csv(primary)
        self.codebook = pd.read_csv(codebook, index_col = 0)
        self.pseudocolors = self.codebook.iloc[:,1].max()
        self.readouts = pd.read_csv(readouts)
    
    def construct_seqfish_probes(self, progress_callback=None):

        #isolate primary that matched codebook
        primary = self.primary.loc[self.primary['<GENE>'].isin(self.codebook.index)].reset_index(drop=True)
        #only take gene and probe sequence
        primary = primary[["<GENE>", "<PROBE>"]]
        #list of channels
        channel_names = ["647","561","488"]
        #isolate barcodes
        barcodes = []
        #total pseudocolor / number of channels will be hybs per round
        hyb_per_round = int(self.pseudocolors/3)
        num_tasks = len(primary)
        for i in range(len(primary)):
            if progress_callback is not None:
                progress = 100 * (i + 1) / num_tasks
                progress_callback(progress)
            #get gene name from probes csv
            genename = primary.iloc[i]["<GENE>"]
            #get pseudocolor info from codebook
            pseudocolors = self.codebook[self.codebook.index == genename].values[0]
            #get the amplifiers that match pseudocolors
            amplifiers_for_gene = []
            round_info = 0
            for ps in pseudocolors:
                channel = (ps%3)-1
                if channel == -1:
                    channel = 3-1
                #isolate desired channel amplifiers
                ch_amplifiers = self.readouts[self.readouts.channels.str.contains(channel_names[channel])]["amplifier_sequence"].values
                #counts the spacing to get hyb info
                hyb = len(np.arange(channel+1,ps+1,3))-1
                amplifiers_for_gene.append(ch_amplifiers[hyb+round_info])
                round_info += hyb_per_round
            barcodes.append(amplifiers_for_gene)
        
        #convert to df
        final_codes = pd.DataFrame(barcodes)  
        final_codes.columns = [f"Readout_{i+1}" for i in range(len(self.codebook.columns))]  
        #add in barcode sequence info
        primary = pd.concat([primary, final_codes],axis=1)   
        #add in fwd and reverse primer sequences along with bridge sequence
        primary["Fwd_Primer"] = "GCCCCATCATGTGCCTTTCc"
        primary["Rev_Primer"] = "cTATAGTGAGTCGTATTACCGGCC"
        primary["Bridge"] = "CACTCTACG"
        #add spacers between readout sites
        for i in range(len(self.codebook.columns)-2):
            primary[f"Spacer{i+1}"] = "ctatac"
        #add spacers between primary
        for i in range(2):
            primary[f"Primary_Spacer{i+1}"] = "caac"
        #grab indicies for reorganizing
        columns = primary.columns
        fwd_idx = np.argwhere(columns == "Fwd_Primer")[0][0]
        rev_idx = np.argwhere(columns == "Rev_Primer")[0][0]
        primary_idx = np.argwhere(columns == "<PROBE>")[0][0]
        readout_idxs = [np.argwhere(columns == f"Readout_{i+1}")[0][0] for i in range(len(self.codebook.columns))]
        spacer_idxs = [np.argwhere(columns == f"Spacer{i+1}")[0][0] for i in range(len(self.codebook.columns)-2)]
        primary_spacer_idxs = [np.argwhere(columns == f"Primary_Spacer{i+1}")[0][0] for i in range(2)]
        bridge_idx = np.argwhere(columns == "Bridge")[0][0]
        #create new index array
        final_indicies = []
        add_primary = int(len(readout_idxs)/2)
        spacer_num = 0
        for i in range(len(readout_idxs)):
            if i == 0:
                #add forward
                final_indicies.append(fwd_idx)
            if i == len(readout_idxs)-1:
                #add final readout
                final_indicies.append(readout_idxs[i])
                #add bridge site
                final_indicies.append(bridge_idx)
                #add reverse
                final_indicies.append(rev_idx)
                break
            if i != add_primary:
                #add readout
                final_indicies.append(readout_idxs[i])
                #add readout spacer
                final_indicies.append(spacer_idxs[spacer_num])
                spacer_num += 1
            else:
                #remove last index
                final_indicies.pop()
                spacer_num -= 1
                #add primary spacer
                final_indicies.append(primary_spacer_idxs[0])
                #add primary
                final_indicies.append(primary_idx)
                #add primary spacer
                final_indicies.append(primary_spacer_idxs[1])
                #add readout
                final_indicies.append(readout_idxs[i])
                #add readout spacer
                final_indicies.append(spacer_idxs[spacer_num])
                spacer_num += 1
            
        #reshuffle columns
        primary = primary.set_index("<GENE>")
        primary = primary[columns[final_indicies]]
        #assemble probes
        complete_probe = ["".join(sublist) for sublist in primary.values.tolist()]
        primary["Final_Probe"] = complete_probe

        #write out final probes
        output_dir = Path.home() / Path("ProbeBuilder/output")
        output_dir.mkdir(parents=True, exist_ok = True)
        primary.to_csv(str(output_dir / "assembled_inverted_padlocks.csv"))
        #write out primers and bridges required
        with open(str(output_dir / "primers_splints_to_buy.txt"), "w+") as f:
            f.write("Order the following sequences if you do not have it already.\n")
            f.write("Forward Primer: GCCCCATCATGTGCCTTTCC \n")
            f.write("Reverse Primer: GGCCGGTAATACGACTCACTATAG \n")
            f.write("Reverse Transcription Primer: /5Phos/GCCCCATCATGTGCCTTTCC \n")
            f.write("LNA Splint: C+AT+GAT GGG GCG CGT AG+AGT+G \n")
            f.close()

    