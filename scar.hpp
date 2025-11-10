#ifndef SCAR_HPP
#define SCAR_HPP

#include <htslib/sam.h>
#include <htslib/hts.h>
#include <string>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <memory>

struct MetadataRecord {
    std::string barcode;
    std::string pool_id;
    std::string donor_id;
    std::string celltype_lev1;
    std::string celltype_lev2;
    std::string bamfile;
    std::string bed_file;
};

struct CellTypeDonorKey {
    std::string celltype;
    std::string donor_id;
    
    bool operator<(const CellTypeDonorKey& other) const {
        if (celltype != other.celltype) return celltype < other.celltype;
        return donor_id < other.donor_id;
    }
};

struct BamGroupInfo {
    std::string bamfile;
    std::string donor_id;
    std::string celltype;
    std::string bed_file;
    std::unordered_set<std::string> barcodes;
    size_t cell_count;
};

struct Interval {
    int32_t start;
    int32_t end;
};

class BedRegions {
public:
    bool loadBed(const std::string& filename);
    bool overlaps(const std::string& chrom, int32_t pos, int32_t len) const;
    size_t getIntervalCount() const;
    
private:
    std::map<std::string, std::vector<Interval>> region_intervals_;
};

struct OutputFile {
    std::string filename;
    samFile* handle;
    std::vector<bam1_t*> write_buffer;
    std::shared_ptr<BedRegions> bed_regions;
    size_t written_count;
    
    static constexpr size_t BUFFER_SIZE = 10000;
    
    bool flush(sam_hdr_t* header) {
        for (auto* b : write_buffer) {
            if (sam_write1(handle, header, b) < 0) {
                return false;
            }
            bam_destroy1(b);
        }
        write_buffer.clear();
        return true;
    }
};

class MetadataProcessor {
public:
    MetadataProcessor(int min_cells = 10);
    
    bool loadMetadata(const std::string& filename, const std::string& celltype_column);
    bool loadPeakMatrix(const std::string& peak_matrix_file, const std::string& output_bed_dir);
    void filterMetadata();
    bool writeFilteredMetadata(const std::string& output_filename);
    
    const std::map<CellTypeDonorKey, BamGroupInfo>& getFilteredGroups() const {
        return filtered_groups_;
    }
    
private:
    int min_cells_;
    std::string celltype_column_;
    std::string bed_output_dir_;
    std::vector<MetadataRecord> records_;
    std::map<CellTypeDonorKey, BamGroupInfo> filtered_groups_;
    std::map<std::string, std::string> celltype_to_bed_;
    
    std::string getCellType(const MetadataRecord& record) const;
    void resolvePoolConflicts();
};

class BamRouter {
public:
    BamRouter();
    ~BamRouter();
    
    void setRemoveDuplicates(bool remove) { remove_duplicates_ = remove; }
    
    bool routeBamSinglePass(const std::string& input_bam,
                           const std::map<CellTypeDonorKey, BamGroupInfo>& groups,
                           const std::string& output_dir);
    
private:
    bool remove_duplicates_;
};

#endif // SCAR_HPP

