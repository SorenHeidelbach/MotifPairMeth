
use ahash::{HashMap, HashMapExt};
use motif_methylation_state::utils::modtype::ModType;
use motif_methylation_state::utils::strand::Strand;
use crate::pileup::{PileupRecord, PileupChunk};
use crate::motif::Motif;
use regex::Regex;

#[derive(Debug, Clone)]
pub struct Contig {
    pub reference: String,
    pub sequence: String,
    pub records: HashMap<(usize, Strand, ModType), PileupRecord>,
}

impl Contig {
    pub fn new(reference: &str, sequence: &str) -> Self {
        Self {
            reference: reference.to_string(),
            sequence: sequence.to_string(),
            records: HashMap::new(),
        }
    }

    pub fn add_record(&mut self, record: PileupRecord) {
        let key = (record.position, record.strand, record.mod_type);
        self.records.insert(key, record);
    }

    pub fn add_records(&mut self, records: PileupChunk) {
        if records.reference != self.reference {
            panic!("Reference mismatch: {} != {}", records.reference, self.reference);
        }
        self.records.reserve(records.records.len());
        for record in records.records {
            self.add_record(record);
        }
    }

    pub fn find_motif_indeces(&self, motif: Motif) -> Vec<usize> {
        let mut indices = Vec::new();
        let motif_regex = motif.regex().unwrap();
        let re = Regex::new(&motif_regex).unwrap();
        re.find_iter(&self.sequence)
            .map(|m| indices.push(m.start() as usize + motif.position as usize));
        indices
    }

    pub fn find_complement_motif_indeces(&self, motif: Motif) -> Vec<usize> {
        let mut indices = Vec::new();
        let complement_motif = motif.reverse_complement().unwrap();
        let motif_regex = complement_motif.regex().unwrap();
        let re = Regex::new(&motif_regex).unwrap();
        re.find_iter(&self.sequence)
            .map(|m| indices.push(m.start() as usize + complement_motif.position as usize));
        indices
    }

    pub fn get_records(&self, position: Vec<usize>, strand: Strand, mod_type: ModType) ->  Vec<&PileupRecord> {
        position
            .iter()
            .map(|&p| {
                let key = (p, strand, mod_type);
                self.records.get(&key).unwrap()
            })
            .collect()
    }
}





