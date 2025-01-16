use crate::motif::{Motif, MotifLike};
use crate::pileup::{PileupChunk, PileupRecord};
use ahash::{HashMap, HashMapExt};
use memopair::utils::modtype::ModType;
use memopair::utils::strand::Strand;
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
            panic!(
                "Reference mismatch: {} != {}",
                records.reference, self.reference
            );
        }
        self.records.reserve(records.records.len());
        for record in records.records {
            self.add_record(record);
        }
    }

    pub fn find_motif_indeces(&self, motif: &Motif) -> Option<Vec<usize>> {
        let mut indices = Vec::new();
        let motif_regex = motif.regex().unwrap();
        let re = Regex::new(&motif_regex).unwrap();
        // Find matches in the contig sequence of the motif
        re.find_iter(&self.sequence)
            .map(|m| indices.push(m.start() as usize + motif.position as usize))
            .for_each(drop);
        if indices.is_empty() {
            return None;
        }
        Some(indices)
    }

    pub fn find_complement_motif_indeces(&self, motif: &Motif) -> Option<Vec<usize>> {
        let mut indices = Vec::new();
        let complement_motif = motif.reverse_complement().unwrap();
        let motif_regex = complement_motif.regex().unwrap();
        let re = Regex::new(&motif_regex).unwrap();
        re.find_iter(&self.sequence)
            .map(|m| indices.push(m.start() as usize + complement_motif.position as usize))
            .for_each(drop);
        if indices.is_empty() {
            return None;
        }
        Some(indices)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::motif::Motif;
    use crate::pileup::PileupRecord;
    use memopair::utils::modtype::ModType;
    use memopair::utils::strand::Strand;

    #[test]
    fn test_contig_add_record() {
        let mut contig = Contig::new("test", "ACGT");
        let record = PileupRecord {
            reference: "test".to_string(),
            position: 0,
            strand: Strand::Positive,
            mod_type: ModType::SixMA,
            n_mod: 1,
            n_valid_cov: 1,
            n_canonical: 1,
            n_diff: 0,
        };
        contig.add_record(record.clone());
        assert_eq!(contig.records.len(), 1);
        assert_eq!(
            contig.records.get(&(0, Strand::Positive, ModType::SixMA)),
            Some(&record)
        );
    }

    #[test]
    fn test_contig_add_records() {
        let mut contig = Contig::new("test", "ACGT");
        let records = PileupChunk {
            reference: "test".to_string(),
            records: vec![
                PileupRecord {
                    reference: "test".to_string(),
                    position: 0,
                    strand: Strand::Positive,
                    mod_type: ModType::SixMA,
                    n_mod: 1,
                    n_valid_cov: 1,
                    n_canonical: 1,
                    n_diff: 0,
                },
                PileupRecord {
                    reference: "test".to_string(),
                    position: 1,
                    strand: Strand::Positive,
                    mod_type: ModType::SixMA,
                    n_mod: 1,
                    n_valid_cov: 1,
                    n_canonical: 1,
                    n_diff: 0,
                },
            ],
        };
        contig.add_records(records.clone());
        assert_eq!(contig.records.len(), 2);
        assert_eq!(
            contig.records.get(&(0, Strand::Positive, ModType::SixMA)),
            Some(&records.records[0])
        );
        assert_eq!(
            contig.records.get(&(1, Strand::Positive, ModType::SixMA)),
            Some(&records.records[1])
        );
    }

    #[test]
    fn test_contig_find_motif_indeces() {
        let contig = Contig::new("test", "ACGTACGTACGTACGT");
        let motif = Motif::new("ACGT", "6mA", 0).unwrap();
        let indeces = contig.find_motif_indeces(&motif).unwrap();
        assert_eq!(indeces, vec![0, 4, 8, 12]);

        let contig = Contig::new(
            "test",
            "ACCCCGGAGGTCGTACGCCGGATCCGGTACCGGACGTACCGGTCGCCGGAT",
        );
        let motif = Motif::new("CCGGA", "6mA", 4).unwrap();
        let indeces = contig.find_motif_indeces(&motif);
        assert_eq!(indeces, Some(vec![7, 21, 33, 49]));
    }

    #[test]
    fn test_contig_find_complement_motif_indeces() {
        let contig = Contig::new("test", "ACGTACGTACGTACGT");
        let motif = Motif::new("ACGT", "6mA", 0).unwrap();
        let indeces = contig.find_complement_motif_indeces(&motif);
        assert_eq!(indeces, Some(vec![3, 7, 11, 15]));

        let contig = Contig::new(
            "test",
            "ACCTCCGGCCGGAGGTCGTACGCCGGATCCGGTCCGGTCCGGTACCGGACGTACCGGTCGCCGGAT",
        );
        let motif = Motif::new("CCGGA", "6mA", 4).unwrap();
        let indeces = contig.find_complement_motif_indeces(&motif);
        assert_eq!(indeces, Some(vec![3, 27, 32, 37]));

        let contig = Contig::new("test", "CCTCCTCCTCCTCCTCC");
        let motif = Motif::new("CCTCC", "5mC", 0).unwrap();
        let indeces = contig.find_complement_motif_indeces(&motif);
        assert_eq!(indeces, None);

        let contig = Contig::new("test", "GGAGGAGGAGGAGGAGG");
        let motif = Motif::new("CCTCC", "5mC", 0).unwrap();
        let indeces = contig.find_complement_motif_indeces(&motif);
        assert_eq!(indeces, Some(vec![4, 10, 16])); // only count full matches
    }
}
