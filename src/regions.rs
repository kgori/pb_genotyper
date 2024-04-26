use std::collections::HashMap;

#[derive(Debug, Clone)]
pub struct Region {
    pub seqname: String,
    pub start: i64,
    pub end: i64,
}

pub fn sort_regions_in_place(regions: &mut [Region]) {
    regions.sort_by(|a, b| {
        a.seqname
            .cmp(&b.seqname)
            .then_with(|| a.start.cmp(&b.start))
            .then_with(|| a.end.cmp(&b.end))
    });
}

pub fn compress_regions(regions: &[Region]) -> Vec<Region> {
    let mut compressed = vec![];
    let mut current = regions[0].clone();
    for region in regions.iter().skip(1) {
        if region.seqname == current.seqname && region.start <= current.end {
            current.end = region.end;
        } else {
            compressed.push(current.clone());
            current = region.clone();
        }
    }
    compressed.push(current);
    compressed
}

// Converts a vector of regions into a hashmap, where the key is the
// chromosome name and the value is a sorted vector of regions on that chromosome.
pub fn convert_regions_vec_to_hashmap(regions: Vec<Region>) -> HashMap<String, Vec<Region>> {
    let mut regions_map = HashMap::new();
    for region in regions {
        let entry = regions_map.entry(region.seqname.clone()).or_insert(vec![]);
        entry.push(region);
    }
    for (_, regions) in regions_map.iter_mut() {
        regions.sort_by(|a, b| a.start.cmp(&b.start).then_with(|| a.end.cmp(&b.end)));
    }
    regions_map
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sort_regions_in_place_different_seqnames() {
        let mut regions = vec![
            Region {
                seqname: "chr2".to_string(),
                start: 100,
                end: 200,
            },
            Region {
                seqname: "chr1".to_string(),
                start: 100,
                end: 200,
            },
        ];
        sort_regions_in_place(&mut regions);
        assert_eq!(regions[0].seqname, "chr1");
        assert_eq!(regions[1].seqname, "chr2");
    }

    #[test]
    fn test_sort_regions_in_place_same_seqname_different_starts() {
        let mut regions = vec![
            Region {
                seqname: "chr1".to_string(),
                start: 200,
                end: 300,
            },
            Region {
                seqname: "chr1".to_string(),
                start: 100,
                end: 200,
            },
        ];
        sort_regions_in_place(&mut regions);
        assert_eq!(regions[0].start, 100);
        assert_eq!(regions[1].start, 200);
    }

    #[test]
    fn test_sort_regions_in_place_same_seqname_same_starts() {
        let mut regions = vec![
            Region {
                seqname: "chr1".to_string(),
                start: 100,
                end: 300,
            },
            Region {
                seqname: "chr1".to_string(),
                start: 100,
                end: 200,
            },
        ];
        sort_regions_in_place(&mut regions);
        assert_eq!(regions[0].end, 200);
        assert_eq!(regions[1].end, 300);
    }

    #[test]
    fn test_compress_regions_non_overlapping() {
        let regions = vec![
            Region {
                seqname: "chr1".to_string(),
                start: 100,
                end: 200,
            },
            Region {
                seqname: "chr1".to_string(),
                start: 300,
                end: 400,
            },
        ];
        let compressed = compress_regions(&regions);
        assert_eq!(compressed.len(), 2);
        assert_eq!(compressed[0].start, 100);
        assert_eq!(compressed[1].start, 300);
    }

    #[test]
    fn test_compress_regions_overlapping() {
        let regions = vec![
            Region {
                seqname: "chr1".to_string(),
                start: 100,
                end: 200,
            },
            Region {
                seqname: "chr1".to_string(),
                start: 150,
                end: 250,
            },
            Region {
                seqname: "chr1".to_string(),
                start: 240,
                end: 300,
            },
        ];
        let compressed = compress_regions(&regions);
        assert_eq!(compressed.len(), 1);
        assert_eq!(compressed[0].start, 100);
        assert_eq!(compressed[0].end, 300);
    }

    #[test]
    fn test_compress_regions_multiple_seqnames() {
        let regions = vec![
            Region {
                seqname: "chr1".to_string(),
                start: 100,
                end: 200,
            },
            Region {
                seqname: "chr2".to_string(),
                start: 100,
                end: 200,
            },
        ];
        let compressed = compress_regions(&regions);
        assert_eq!(compressed.len(), 2);
        assert_eq!(compressed[0].seqname, "chr1");
        assert_eq!(compressed[1].seqname, "chr2");
    }

    #[test]
    fn test_convert_regions_vec_to_hashmap() {
        let regions = vec![
            Region {
                seqname: "chr1".to_string(),
                start: 100,
                end: 200,
            },
            Region {
                seqname: "chr1".to_string(),
                start: 150,
                end: 250,
            },
            Region {
                seqname: "chr2".to_string(),
                start: 100,
                end: 200,
            },
        ];
        let regions_map = convert_regions_vec_to_hashmap(regions);
        assert_eq!(regions_map.len(), 2);
        assert_eq!(regions_map.get("chr1").unwrap().len(), 2);
        assert_eq!(regions_map.get("chr2").unwrap().len(), 1);
    }
}
