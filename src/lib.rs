#![allow(unused, non_snake_case)]


use hashbrown::HashMap;
use itertools::Combinations;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use serde::{Deserialize, Serialize};
use serde_json::from_str;
use std::fs::read_to_string;
use std::ops::DerefMut;
use std::{path::PathBuf, thread, panic};
use std::sync::atomic::{AtomicUsize, AtomicBool, Ordering};
use std::time::Duration;
use std::process;

static GLOBAL_DO_PANIC: AtomicBool = AtomicBool::new(false);



fn busy_wait_panic_thread() {
    thread::spawn(|| {
        loop {
            match GLOBAL_DO_PANIC.load(Ordering::Relaxed) {
                true => {
                    println!("Exiting prematurely...");
                    process::abort();
                },
                false => (),
            }
        }
    });
}

fn panic_silently(message: &str) {
    println!("\nERROR CAUGHT: {}\n", message);

    GLOBAL_DO_PANIC.store(true, Ordering::Relaxed);   
}


#[derive(Serialize, Deserialize, Clone, Debug)]
struct GeneTable {
    table_id: String,
    A: Option<Vec<String>>,
    L: Option<Vec<String>>,
    W: Option<Vec<String>>,
    Q: Option<Vec<String>>,
    Y: Option<Vec<String>>,
    E: Option<Vec<String>>,
    C: Option<Vec<String>>,
    D: Option<Vec<String>>,
    F: Option<Vec<String>>,
    G: Option<Vec<String>>,
    H: Option<Vec<String>>,
    I: Option<Vec<String>>,
    M: Option<Vec<String>>,
    K: Option<Vec<String>>,
    P: Option<Vec<String>>,
    R: Option<Vec<String>>,
    S: Option<Vec<String>>,
    V: Option<Vec<String>>,
    N: Option<Vec<String>>,
    T: Option<Vec<String>>,
    X: Option<Vec<String>>,
    B: Option<Vec<String>>,
    Z: Option<Vec<String>>,
    J: Option<Vec<String>>,
    STOP: Option<Vec<String>>,
    DASH: Option<Vec<String>>,
}

impl GeneTable {
    pub fn from(s: String, table_id: u8) -> Self {
        let gene_tables: Vec<GeneTable> = from_str(&s).unwrap();

        let mut gene_table = gene_tables
            .iter()
            .cloned()
            .filter(|x| x.table_id == table_id.to_string())
            .collect::<Vec<GeneTable>>()
            .get(0)
            .unwrap()
            .clone();

        gene_table.set_special();

        gene_table
    }

    pub fn set_special(&mut self) {
        self.DASH = Some(vec!["---".to_string()]);
        self.X = Some(vec![
            'A', 
            'L',
            'W',
            'Q',
            'Y',
            'E',
            'C',
            'D',
            'F',
            'G',
            'H',
            'I',
            'M',
            'K',
            'P',
            'R',
            'S',
            'V',
            'N',
            'T',
            '*',
        ].iter().cloned()
                .map(|c| {
                    match self.get(&c) {
                        Some(t) => t.clone(),
                        None => vec![],
                    }
                })
                .chain(vec![vec![
                        "NNN".to_string(), 
                        "NNG".to_string(), 
                        "NNA".to_string(), 
                        "NAC".to_string(), 
                        "NAC".to_string(), 
                        "ANC".to_string(), 
                        "ANN".to_string()
                    ]]                    
                    .into_iter()
                )
                .flatten()
                .collect::<Vec<String>>()
            );
        self.B = Some(vec![self.D.clone(), self.N.clone()]
                    .into_iter()
                    .map(|v| v.into_iter().flatten())
                    .flatten()
                    .collect::<Vec<String>>()
                );

        self.Z = Some(vec![self.E.clone(), self.Q.clone()]
                    .into_iter()
                    .map(|v| v.into_iter().flatten())
                    .flatten()
                    .collect::<Vec<String>>()
            );
        self.J = Some(vec![self.I.clone(), self.L.clone()]
                    .into_iter()
                    .map(|v| v.into_iter().flatten())
                    .flatten()
                    .collect::<Vec<String>>()
        );
            
    }

    pub fn get(&self, c: &char) -> &Option<Vec<String>> {
        match c {
            'A' => &self.A,
            'L' => &self.L,
            'W' => &self.W,
            'Q' => &self.Q,
            'Y' => &self.Y,
            'E' => &self.E,
            'C' => &self.C,
            'D' => &self.D,
            'F' => &self.F,
            'G' => &self.G,
            'H' => &self.H,
            'I' => &self.I,
            'M' => &self.M,
            'K' => &self.K,
            'P' => &self.P,
            'R' => &self.R,
            'S' => &self.S,
            'V' => &self.V,
            'N' => &self.N,
            'T' => &self.T,
            '*' => &self.STOP,
            '-' => &self.DASH,
            'X' => &self.X,
            'B' => &self.B,
            'J' => &self.J,
            'Z' => &self.Z,
            _ => {
                panic_silently("No acceptable index for gene table, must be one of the 18 AA chars + 4 ambiguous chars.");
                &self.STOP 
            }
        }
    }

    }

impl std::ops::Index<char> for GeneTable {
    type Output = Option<Vec<String>>;

    fn index(&self, index: char) -> &Self::Output {
        self.get(&index)
    }
}



struct AminoAcidTranslator((String, String), (String, String));

impl AminoAcidTranslator {
    pub fn do_checks(&mut self) {
        let AminoAcidTranslator((aa_header, aa), (nt_header, nt)) = self;

        if aa_header != nt_header {
            panic_silently(&format!("AA header -> {} is not the same as NT header -> {}", aa_header, nt_header));
        }

        let len_aa = aa.len();
        let len_nt = nt.len();
        let aa_filt_mul = aa.chars().filter(|c| *c != '-').count() * 3;

        if len_nt != aa_filt_mul {
            let longer_shorter = match aa_filt_mul > len_nt {
                true => (
                    format!("(AA -> {})", aa_header),
                    format!("(NT -> {})", nt_header)
                ),
                false => (
                    format!("(NT -> {})", nt_header),
                    format!("(AA -> {})", aa_header),
                )
            };

            let diff = {
                let num_marker = match aa_filt_mul > len_nt {
                        true => ((aa_filt_mul - len_nt) / 3, "PEP char(s)"),
                        false => ((len_nt - aa_filt_mul) / 3, "NT triplet(s)"),
                    };
                format!("with a difference of {} {}", num_marker.0, num_marker.1)
            };

            panic_silently(&format!("{} is larger than {} {}", longer_shorter.0, longer_shorter.1, diff));
            
        }                     
    }

    
    pub fn streamline_amino_acid(&mut self) {
        let AminoAcidTranslator((header, amino_acid), _) = self;

        let mut amino_acid_trimmed = amino_acid.trim().to_uppercase();

        let mut amino_acid_filtered = String::new();

        amino_acid_trimmed
            .char_indices()
            .for_each(|(i, c)| {
                match !vec![
                    'A', 
                    'L',
                    'W',
                    'Q',
                    'Y',
                    'E',
                    'C',
                    'D',
                    'F',
                    'G',
                    'H',
                    'I',
                    'M',
                    'K',
                    'P',
                    'R',
                    'S',
                    'V',
                    'N',
                    'T',
                    '*',
                    '-',
                ].contains(&c) {
                    true => {                       
                        amino_acid_filtered.push('X');
                    },
                    false => amino_acid_filtered.push(c),
                } 
        });

        *amino_acid = amino_acid_filtered;

    }

    fn error_out(&self) {
        let AminoAcidTranslator((header, amino_acid), (_, compare_dna)) = self;


        panic_silently(&format!(r#" 
======
MISMATCH ERROR:
    The following Amino Acid failed to match with its source Nucleotide pair.

        Header: `{}`,                    
        ===
        Amino Acid: `{}`,
        ===
        Source Nucleotide: `{}`,
=======
            "#, 
                header, 
                amino_acid, 
                compare_dna,
            ));
        }   

    pub fn reverse_translate_and_compare(
        &self,
        gene_table: &GeneTable,
    ) -> String {
        let AminoAcidTranslator((header, amino_acid), (_, compare_dna)) = self;

        let mut compare_triplets = (0..compare_dna.len())
                                                .step_by(3)
                                                .map(|i| compare_dna[i..i + 3].to_string())
                                                .into_iter();

        let final_taxon = amino_acid
            .chars()
            .map(|aa| {
                match aa == '-' {
                    true => "---".to_string(),
                    false => {
                        match aa.is_ascii_digit() {
                            true => {
                                let d = aa.to_digit(110).unwrap();

                                ".".repeat(d as usize).to_string()
                            },
                            false => {
                                let mut taxa_triplets = gene_table.get(&aa);

                                match taxa_triplets {
                                    Some(taxa) => {
                                        let mut taxa_mut = taxa.clone();
                                        let original_triplet = compare_triplets.next().unwrap();

                                        taxa_mut.retain(|s| s == &original_triplet);
                                        
                                        match taxa_mut.get(0) {
                                            Some(t) => t.clone(),
                                            None => {
                                                //self.error_out();
                                                "".to_string()
                                            },
                                        }
                                    },
                                    None => {
                                        panic_silently("Genetic table does not have the pep. Perhaps you've chosen the wrong table index?");
                                        "".to_string()
                                        
                                    }
                                }
                            }
                        }
                    }
                }

                
            })
            .collect::<Vec<String>>()
            .join("");     
           
        final_taxon
    }
}

#[pyclass]
struct PEP2NT(GeneTable);

#[pymethods]
impl PEP2NT {
    #[new]
    pub fn new(s: String, table_id: u8) -> Self {
        busy_wait_panic_thread();
    
        let gene_table = GeneTable::from(s, table_id);

        Self(gene_table)
    }   

    pub fn pro2codon( 
        &self,
        amino_seqs: Vec<(String, String)>,
        nuc_seqs: Vec<(String, String)>,
    ) -> Vec<(String, String)>
    {
        let PEP2NT(gene_table) = self;

        let aa_seq_len = amino_seqs.len();
        let nt_seq_len = nuc_seqs.len();

        if aa_seq_len != nt_seq_len {
            let longer_shorter = match aa_seq_len > nt_seq_len {
                true => ("AA", "NT"),
                false => ("NT", "AA")
            };

            let diff = match aa_seq_len > nt_seq_len {
                true => (aa_seq_len as isize - aa_seq_len as isize).abs(),
                false => (nt_seq_len as isize - aa_seq_len as isize).abs(),
            };

            panic_silently(&format!("Length of the {} sequence is longer than the length of {} sequence by a number of {}.", longer_shorter.0, longer_shorter.1, diff));
        
        }
        

        amino_seqs
            .iter()
            .cloned()
            .zip(nuc_seqs.iter().cloned())
            .map(|((aa_header, aa), (nt_header, nt))| {
                if GLOBAL_DO_PANIC.load(Ordering::Relaxed) {
                    thread::sleep(Duration::from_secs(5));
                }

                let mut amino_acid = AminoAcidTranslator(
                    (aa_header.clone(), aa.clone()),
                    (nt_header.clone(), nt.clone()),
                );
                amino_acid.do_checks();
                amino_acid.streamline_amino_acid();

                let codon = amino_acid.reverse_translate_and_compare(&gene_table);
                
                (nt_header.clone(), codon)
            })
            .collect()          
    }
}



/// A Python module implemented in Rust.
#[pymodule]
fn pro2codon(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PEP2NT>();
    Ok(())
}


#[test]
fn test_revtrans_success() {
    let s = r#"
    [{
        "table_id": "1",
        "F": [
            "TTT",
            "TTC"
        ],
        "L": [
            "TTA",
            "TTG",
            "CTT",
            "CTC",
            "CTA",
            "CTG"
        ],
        "S": [
            "TCT",
            "TCC",
            "TCA",
            "TCG",
            "AGT",
            "AGC"
        ],
        "Y": [
            "TAT",
            "TAC"
        ],
        "*": [
            "TAA",
            "TAG",
            "TGA"
        ],
        "C": [
            "TGT",
            "TGC"
        ],
        "W": [
            "TGG"
        ],
        "P": [
            "CCT",
            "CCC",
            "CCA",
            "CCG"
        ],
        "H": [
            "CAT",
            "CAC"
        ],
        "Q": [
            "CAA",
            "CAG"
        ],
        "R": [
            "CGT",
            "CGC",
            "CGA",
            "CGG",
            "AGA",
            "AGG"
        ],
        "I": [
            "ATT",
            "ATC",
            "ATA"
        ],
        "M": [
            "ATG"
        ],
        "T": [
            "ACT",
            "ACC",
            "ACA",
            "ACG"
        ],
        "N": [
            "AAT",
            "AAC"
        ],
        "K": [
            "AAA",
            "AAG"
        ],
        "V": [
            "GTT",
            "GTC",
            "GTA",
            "GTG"
        ],
        "A": [
            "GCT",
            "GCC",
            "GCA",
            "GCG"
        ],
        "D": [
            "GAT",
            "GAC"
        ],
        "E": [
            "GAA",
            "GAG"
        ],
        "G": [
            "GGT",
            "GGC",
            "GGA",
            "GGG"
        ]
    }
    
]
    
    "#.to_string();
    let table_id = 1u8;

    let pep2nt = PEP2NT::new(s, table_id);

    let nuc = vec![("aa".to_string(), "ATGTGCCATCAGCATCCGAATCATTTGCTGAAGCAGTCTCTCCAAACGTACGTGTTCCGCCGAAATGAGCTACTGGACGCCGGCAAAAGTGCACTGCGGCATCGCATTCGCCTTCAGGAGCGAAATTGGATGACGACAGGGGCCAAACCACAGCAGCAGCAAGTCGTGCGAACCGGGGAAAGGAGCTCAAGGGTTCTGAGGTTACTGCTGGCAGGAACCGCCGTGAGTGTCACCGTAGGGTATCAATGCCGCCGGTGGTTCGTGCATTGTGAAGCGGCTGTTGTAAACAACCGGCTGGTGGGACAGAAGGTTCCGGGAGCAGGAGATGGGATCAAATTTGACTGGAGGAAATTTTGGAGCTATTTGAGACCGCACCTGGTGAAACTGATTGGGGCGATTATGGCGGCGCTGGCCGTGGCCTATTTTAATATTCAGATTCCCAACATGCTGGGAGTAGTGGTGAATACGCTGTCGAAATATGCCAGGACGAGTTTGAGCGACATTGATTCTTCAACGTTCATGAACGAAATGAAACTGCCATCGCTGCGACTGTTCGGTATGTACATTGCCCAAGCCGGATTCACATTCGTCTACATTTTGCTGCTGAGTCAGATAGGCGAACAGATGGCCGCGAAGATCCGACAGGATCTGTTCAAGCAGATCATCATCCAGGATCTGGAGTTCTTCGACGAAAACCGAACCGGAGAGTTGGTCAACCGGTTGACGGCCGACGTGCAAGACTTCAAATCTAGCTTCAAGCAATGCATCTCGCAAGGATTACGGTCATTTGCCCAGCTGATCGGAGGTGGCATATCGTTGTTCCTTATCTCGCCACAGCTGGCCAGTATCGCTTTGGTATCGGTCCCTGCCGCGGTGGCAATGTTCTCCTTTCTTGGAAAGTCTTTGAGAGCGCTGAGCAAGAAGAGTCAAGCCCAATCGGAACGAGCGACTTCAGTGAGTGAGGAAGCGCTGTCCAACATCAGAACTGTCCGCTCCAGTGCATGCGAGTTTGCCGAAGTTGAACTGTTGCGGCAGGAAACGGAAAAGGCAGCTGAGCTGTCCCAGCAGCTGGGTGTGAAATCGCACTCTATCATGCCACTGACTAATTTATACGTGAACGGCATGGTGCTAAACGACGCTGGTTCTCGTTCAATAAGTGCCGGAGATTTGATGGCATTCTTGGTTGCATCTCAAGGTGTTCAACGCTCTCTTGCCCAAGGCTCTATTCTGCTTGGATCCGTAATCCGTGGAATGACAGCCGGCACCCGTGTATTCGAATACCTCTCGGTTCAGCCTAAAGTAGATCTCAAATACGGACAAATCATTCCCGAGTCGGAAATTCGTGGCGAAATCCGATTCGAAAACGTGTCTTTCACGTATCCTTCCCGACCCAATCATCCTGTTCTCAAAAACTTTTCGCTTGTCCTGAAACCTGGACAAACTGTGGCCCTGGTGGGAGCCAGTGGCTCAGGGAAATCGACCATTGCTTCGCTTCTGGAGCGATTTTACGAGCCAACCGGTGGCCGCATCACCATAGATGGTTATGAACTCTCACAATTGTCGCCTTATTGGCTCCGAGGCCACCTGATAGGATTCATCGAACAGCAACCGATACTGTTCGGAACGTCCATCTACGAGAACATACGCTACGGCCGTCCAGAGGCGTCACAAGCGGAAGTCCTGGAAGCGGCCAAGCTGTCTCAGTCGCATCAGTTTGTGAGCAAATTGCCTGAGGGCTACGAGACGCCAGTTGGCGAAAGGGGCATCCAACTGAGTGGTGGCCAACGGCAAAGGATAGCAATAGCTAGAGCTTTGCTCAAACAGCCCTCGGTGTTGATCTTGGACGAGGCTACCAGTGCATTAGATGCTTCCAGTGAGGCCATCGTGCAGAAAGCCCTTGATGCAGCGGTAGTCGACAGAACTACGCTGGTTATTGCCCACCGGCTATCAACCATCAGGAACGCGGATGTGATTGTGGTGCTGGAAAAAGGCCGAATAGTAGAGATTGGCAATCATGATGCGCTGCTGCGCAAGAAGGGCTACTACTTTGAACTGGTCAAACAGCAAGAACGGGAACAACGCGAGGAACAACAGCAAAGGGCCTACGGA".to_string())];
    let aa = vec![("aa".to_string(), "MC-------HQHPNHLLKQSLQTYVFRRNELLDAGKSALR-HRIRLQERNWMTTGAKPQQQQVVRTGERSSRVLR-------LLLAGTA----VSVTVGYQCRRWF-------VHCEAAVVNNRLVG--QKVPG-AGDGIKFDWRKFWSYLRPHLVKLIGAIMAALAVAYFNIQIPNMLGVVVNTLSKYARTSLSDIDSSTFMNEMKLPSLRLFGMYIAQAGFTFVYILLLSQIGEQMAAKIRQDLFKQIIIQDLEFFDENRTGELVNRLTADVQDFKSSFKQCISQGLRSFAQLIGGGISLFLISPQLASIALVSVPAAVAMFSFLGKSLRALSKKSQAQ----------SERATSVSEEALSNIRTVRSSACEFAEVELLRQETEKAAELSQQLGVKSHSIMPLTNLYVNGMVLND-------AGSRSISAGDLMAFLVASQGVQRSLAQGSILLGSVIRGMTAGTRVFEYLSVQPKVDLKYGQIIPESEIRGEIRFENVSFTYPSRPNHPVLKNFSLVLKPGQTVALVGASGSGKSTIASLLERFYEPTGGRITIDGYELSQLSPYWLRGHLIGFIEQQPILFGTSIYENIRYGRPEASQAEVLEAAKLSQSHQFVSKLPEGYETPVGERGIQLSGGQRQRIAIARALLKQPSVLILDEATSALDASSEAIVQKALDAAVVDRTTLVIAHRLSTIRNADVIVVLEKGRIVEIGNHDALLRKKGYYFELVKQQEREQREEQQQRAYG------------".to_string())];

    let codon = pep2nt.pro2codon(aa, nuc);

    let should_be = String::from("ATGTGC---------------------CATCAGCATCCGAATCATTTGCTGAAGCAGTCTCTCCAAACGTACGTGTTCCGCCGAAATGAGCTACTGGACGCCGGCAAAAGTGCACTGCGG---CATCGCATTCGCCTTCAGGAGCGAAATTGGATGACGACAGGGGCCAAACCACAGCAGCAGCAAGTCGTGCGAACCGGGGAAAGGAGCTCAAGGGTTCTGAGG---------------------TTACTGCTGGCAGGAACCGCC------------GTGAGTGTCACCGTAGGGTATCAATGCCGCCGGTGGTTC---------------------GTGCATTGTGAAGCGGCTGTTGTAAACAACCGGCTGGTGGGA------CAGAAGGTTCCGGGA---GCAGGAGATGGGATCAAATTTGACTGGAGGAAATTTTGGAGCTATTTGAGACCGCACCTGGTGAAACTGATTGGGGCGATTATGGCGGCGCTGGCCGTGGCCTATTTTAATATTCAGATTCCCAACATGCTGGGAGTAGTGGTGAATACGCTGTCGAAATATGCCAGGACGAGTTTGAGCGACATTGATTCTTCAACGTTCATGAACGAAATGAAACTGCCATCGCTGCGACTGTTCGGTATGTACATTGCCCAAGCCGGATTCACATTCGTCTACATTTTGCTGCTGAGTCAGATAGGCGAACAGATGGCCGCGAAGATCCGACAGGATCTGTTCAAGCAGATCATCATCCAGGATCTGGAGTTCTTCGACGAAAACCGAACCGGAGAGTTGGTCAACCGGTTGACGGCCGACGTGCAAGACTTCAAATCTAGCTTCAAGCAATGCATCTCGCAAGGATTACGGTCATTTGCCCAGCTGATCGGAGGTGGCATATCGTTGTTCCTTATCTCGCCACAGCTGGCCAGTATCGCTTTGGTATCGGTCCCTGCCGCGGTGGCAATGTTCTCCTTTCTTGGAAAGTCTTTGAGAGCGCTGAGCAAGAAGAGTCAAGCCCAA------------------------------TCGGAACGAGCGACTTCAGTGAGTGAGGAAGCGCTGTCCAACATCAGAACTGTCCGCTCCAGTGCATGCGAGTTTGCCGAAGTTGAACTGTTGCGGCAGGAAACGGAAAAGGCAGCTGAGCTGTCCCAGCAGCTGGGTGTGAAATCGCACTCTATCATGCCACTGACTAATTTATACGTGAACGGCATGGTGCTAAACGAC---------------------GCTGGTTCTCGTTCAATAAGTGCCGGAGATTTGATGGCATTCTTGGTTGCATCTCAAGGTGTTCAACGCTCTCTTGCCCAAGGCTCTATTCTGCTTGGATCCGTAATCCGTGGAATGACAGCCGGCACCCGTGTATTCGAATACCTCTCGGTTCAGCCTAAAGTAGATCTCAAATACGGACAAATCATTCCCGAGTCGGAAATTCGTGGCGAAATCCGATTCGAAAACGTGTCTTTCACGTATCCTTCCCGACCCAATCATCCTGTTCTCAAAAACTTTTCGCTTGTCCTGAAACCTGGACAAACTGTGGCCCTGGTGGGAGCCAGTGGCTCAGGGAAATCGACCATTGCTTCGCTTCTGGAGCGATTTTACGAGCCAACCGGTGGCCGCATCACCATAGATGGTTATGAACTCTCACAATTGTCGCCTTATTGGCTCCGAGGCCACCTGATAGGATTCATCGAACAGCAACCGATACTGTTCGGAACGTCCATCTACGAGAACATACGCTACGGCCGTCCAGAGGCGTCACAAGCGGAAGTCCTGGAAGCGGCCAAGCTGTCTCAGTCGCATCAGTTTGTGAGCAAATTGCCTGAGGGCTACGAGACGCCAGTTGGCGAAAGGGGCATCCAACTGAGTGGTGGCCAACGGCAAAGGATAGCAATAGCTAGAGCTTTGCTCAAACAGCCCTCGGTGTTGATCTTGGACGAGGCTACCAGTGCATTAGATGCTTCCAGTGAGGCCATCGTGCAGAAAGCCCTTGATGCAGCGGTAGTCGACAGAACTACGCTGGTTATTGCCCACCGGCTATCAACCATCAGGAACGCGGATGTGATTGTGGTGCTGGAAAAAGGCCGAATAGTAGAGATTGGCAATCATGATGCGCTGCTGCGCAAGAAGGGCTACTACTTTGAACTGGTCAAACAGCAAGAACGGGAACAACGCGAGGAACAACAGCAAAGGGCCTACGGA------------------------------------");

    assert_eq!(codon[0].1, should_be);

}