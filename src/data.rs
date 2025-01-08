use ahash::{HashMap, HashMapExt};
use crate::pileup::{PileupRecord, PileupChunk};
use crate::sequence::Contig;







pub struct GenomeWorkSpaceBuilder {
    pub contigs: HashMap<String, Contig>,
}

impl GenomeWorkSpaceBuilder {
    pub fn new() -> Self {
        Self {
            contigs: HashMap::new(),
        }
    }

    pub fn add_contig(&mut self, contig: &Contig) {
        self.contigs.insert(contig.reference.clone(), contig.clone());
    }

    pub fn add_contigs(&mut self, contigs: &Vec<Contig>) {
        for contig in contigs {
            self.add_contig(contig);
        }
    }

    pub fn push_records(&mut self, records: PileupChunk) {
        let reference = records.reference.clone();
        let contig = self
            .contigs
            .get_mut(&reference)
            .expect("Could not find contig");
        contig.add_records(records);
    }

    pub fn get_contig_ids(&self) -> Option<Vec<&str>> {
        let ids = self.contigs.keys().map(|f|f.as_str()).collect();
        Some(ids)
    }

    pub fn build(self) -> GenomeWorkspace {
        GenomeWorkspace {
            contigs: self.contigs,
        }
    }
}


pub struct GenomeWorkspace {
    pub contigs: HashMap<String, Contig>,
}

impl GenomeWorkspace {
    pub fn get_contig(&self, reference: &str) -> Option<&Contig> {
        self.contigs.get(reference)
    }
}



