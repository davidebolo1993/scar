#include "scar.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <map>
#include <set>
#include <unordered_map>
#include <iomanip>
#include <chrono>
#include <algorithm>

namespace fs = std::filesystem;

// ============================================================================
// BedRegions Implementation
// ============================================================================

bool BedRegions::loadBed(const std::string& filename) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Error: Cannot open BED file: " << filename << std::endl;
        return false;
    }
    
    std::string line;
    size_t line_num = 0;
    size_t intervals_loaded = 0;
    
    while (std::getline(infile, line)) {
        line_num++;
        if (line.empty() || line[0] == '#') continue;
        
        std::istringstream iss(line);
        std::string chrom;
        int32_t start, end;
        
        if (!(iss >> chrom >> start >> end)) {
            std::cerr << "Warning: Skipping malformed BED line " << line_num << std::endl;
            continue;
        }
        
        region_intervals_[chrom].push_back({start, end});
        intervals_loaded++;
    }
    
    for (auto& [chrom, intervals] : region_intervals_) {
        std::sort(intervals.begin(), intervals.end(), 
            [](const Interval& a, const Interval& b) { return a.start < b.start; });
    }
    
    return true;
}

bool BedRegions::overlaps(const std::string& chrom, int32_t pos, int32_t len) const {
    auto it = region_intervals_.find(chrom);
    if (it == region_intervals_.end()) {
        return false;
    }
    
    const auto& intervals = it->second;
    if (intervals.empty()) {
        return false;
    }
    
    int32_t query_end = pos + len;
    
    auto first_candidate = std::lower_bound(intervals.begin(), intervals.end(), pos,
        [](const Interval& a, int32_t p) { return a.end <= p; });
    
    for (auto interval_it = first_candidate; interval_it != intervals.end(); ++interval_it) {
        if (interval_it->start >= query_end) {
            break;
        }
        
        if (interval_it->end > pos) {
            return true;
        }
    }
    
    return false;
}

size_t BedRegions::getIntervalCount() const {
    size_t count = 0;
    for (const auto& [chrom, intervals] : region_intervals_) {
        count += intervals.size();
    }
    return count;
}

// ============================================================================
// MetadataProcessor Implementation
// ============================================================================

MetadataProcessor::MetadataProcessor(int min_cells) : min_cells_(min_cells) {}

std::string MetadataProcessor::getCellType(const MetadataRecord& record) const {
    if (celltype_column_ == "celltypes_LEV1") {
        return record.celltype_lev1;
    } else if (celltype_column_ == "celltypes_LEV2") {
        return record.celltype_lev2;
    }
    return "";
}

bool MetadataProcessor::loadMetadata(const std::string& filename,
                                   const std::string& celltype_column) {
    celltype_column_ = celltype_column;
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Error: Cannot open metadata file: " << filename << std::endl;
        return false;
    }
    
    std::string line;
    bool is_header = true;
    size_t line_num = 0;
    
    while (std::getline(infile, line)) {
        line_num++;
        if (is_header) {
            is_header = false;
            continue;
        }
        
        std::istringstream iss(line);
        MetadataRecord record;
        
        if (!(iss >> record.barcode >> record.pool_id >> record.donor_id
              >> record.celltype_lev1 >> record.celltype_lev2 >> record.bamfile)) {
            std::cerr << "Warning: Skipping malformed line " << line_num << std::endl;
            continue;
        }
        
        records_.push_back(record);
    }
    
    std::cout << "Loaded " << records_.size() << " records from metadata file" << std::endl;
    return true;
}

bool MetadataProcessor::loadPeakMatrix(const std::string& peak_matrix_file, const std::string& output_bed_dir) {
    bed_output_dir_ = output_bed_dir;
    fs::create_directories(output_bed_dir);
    
    std::ifstream infile(peak_matrix_file);
    if (!infile.is_open()) {
        std::cerr << "Error: Cannot open peak matrix file: " << peak_matrix_file << std::endl;
        return false;
    }
    
    std::string header_line;
    std::getline(infile, header_line);
    std::istringstream header_iss(header_line);
    
    std::vector<std::string> celltypes;
    std::string celltype;
    while (header_iss >> celltype) {
        celltypes.push_back(celltype);
    }
    
    std::cout << "Found " << celltypes.size() << " celltypes in peak matrix" << std::endl;
    
    std::map<std::string, std::ofstream> celltype_beds;
    for (const auto& ct : celltypes) {
        std::string bed_path = output_bed_dir + "/" + ct + "_peaks.bed";
        celltype_beds[ct].open(bed_path);
        if (!celltype_beds[ct].is_open()) {
            std::cerr << "Error: Cannot create BED file: " << bed_path << std::endl;
            return false;
        }
        celltype_to_bed_[ct] = bed_path;
    }
    
    std::string line;
    size_t line_num = 1;
    size_t peaks_written = 0;
    
    while (std::getline(infile, line)) {
        line_num++;
        std::istringstream iss(line);
        
        std::string peak_coord;
        iss >> peak_coord;
        
        size_t colon_pos = peak_coord.find(':');
        size_t dash_pos = peak_coord.find('-');
        if (colon_pos == std::string::npos || dash_pos == std::string::npos) {
            std::cerr << "Warning: Malformed peak coordinate at line " << line_num << ": " << peak_coord << std::endl;
            continue;
        }
        
        std::string chrom = peak_coord.substr(0, colon_pos);
        std::string start_str = peak_coord.substr(colon_pos + 1, dash_pos - colon_pos - 1);
        std::string end_str = peak_coord.substr(dash_pos + 1);
        
        for (size_t i = 0; i < celltypes.size(); i++) {
            int has_peak;
            if (!(iss >> has_peak)) {
                std::cerr << "Warning: Missing value for celltype " << celltypes[i] 
                          << " at line " << line_num << std::endl;
                break;
            }
            
            if (has_peak == 1) {
                celltype_beds[celltypes[i]] << chrom << "\t" << start_str << "\t" << end_str << "\n";
                peaks_written++;
            }
        }
    }
    
    for (auto& [ct, bed_file] : celltype_beds) {
        bed_file.close();
    }
    
    std::cout << "Generated " << celltypes.size() << " celltype-specific BED files" << std::endl;
    std::cout << "Total peaks written: " << peaks_written << std::endl;
    
    return true;
}

void MetadataProcessor::resolvePoolConflicts() {
    std::map<CellTypeDonorKey, std::map<std::string, size_t>> pool_counts;
    
    for (const auto& record : records_) {
        std::string celltype = getCellType(record);
        if (celltype.empty()) continue;
        
        CellTypeDonorKey key{celltype, record.donor_id};
        pool_counts[key][record.pool_id]++;
    }
    
    std::map<CellTypeDonorKey, std::string> selected_pools;
    
    for (const auto& [key, pools] : pool_counts) {
        std::string best_pool;
        size_t max_count = 0;
        
        for (const auto& [pool, count] : pools) {
            if (count > max_count) {
                max_count = count;
                best_pool = pool;
            }
        }
        
        if (pools.size() > 1) {
            std::cout << "INFO: Donor " << key.donor_id << " with celltype "
                      << key.celltype << " appears in " << pools.size()
                      << " pools. Selected pool " << best_pool
                      << " with " << max_count << " cells" << std::endl;
        }
        
        selected_pools[key] = best_pool;
    }
    
    for (const auto& record : records_) {
        std::string celltype = getCellType(record);
        if (celltype.empty()) continue;
        
        CellTypeDonorKey key{celltype, record.donor_id};
        if (selected_pools[key] != record.pool_id) continue;
        
        if (filtered_groups_.find(key) == filtered_groups_.end()) {
            std::string bed_file = "";
            auto bed_it = celltype_to_bed_.find(celltype);
            if (bed_it != celltype_to_bed_.end()) {
                bed_file = bed_it->second;
            }
            
            filtered_groups_[key] = BamGroupInfo{
                record.bamfile, record.donor_id, celltype, bed_file, {}, 0
            };
        }
        
        filtered_groups_[key].barcodes.insert(record.barcode);
        filtered_groups_[key].cell_count++;
    }
}

void MetadataProcessor::filterMetadata() {
    resolvePoolConflicts();
    
    size_t total_groups = filtered_groups_.size();
    size_t removed = 0;
    
    for (auto it = filtered_groups_.begin(); it != filtered_groups_.end();) {
        if (it->second.cell_count < static_cast<size_t>(min_cells_)) {
            std::cout << "FILTER: Removing donor " << it->first.donor_id
                      << " with celltype " << it->first.celltype
                      << " (only " << it->second.cell_count << " cells, min="
                      << min_cells_ << ")" << std::endl;
            it = filtered_groups_.erase(it);
            removed++;
        } else {
            ++it;
        }
    }
    
    std::cout << "\nFiltering summary:" << std::endl;
    std::cout << " Total donor/celltype combinations: " << total_groups << std::endl;
    std::cout << " Filtered out (< " << min_cells_ << " cells): " << removed << std::endl;
    std::cout << " Remaining combinations: " << filtered_groups_.size() << std::endl;
}

bool MetadataProcessor::writeFilteredMetadata(const std::string& output_filename) {
    std::ofstream outfile(output_filename);
    if (!outfile.is_open()) {
        std::cerr << "Error: Cannot write to " << output_filename << std::endl;
        return false;
    }
    
    outfile << "barcode\tpool_id\tdonor_id\t" << celltype_column_ << "\tbamfile\tbed_file\n";
    
    size_t total_written = 0;
    for (const auto& record : records_) {
        std::string celltype = getCellType(record);
        if (celltype.empty()) continue;
        
        CellTypeDonorKey key{celltype, record.donor_id};
        
        if (filtered_groups_.find(key) != filtered_groups_.end()) {
            if (filtered_groups_[key].bamfile == record.bamfile) {
                outfile << record.barcode << "\t"
                        << record.pool_id << "\t"
                        << record.donor_id << "\t"
                        << celltype << "\t"
                        << record.bamfile << "\t"
                        << filtered_groups_[key].bed_file << "\n";
                total_written++;
            }
        }
    }
    
    std::cout << "Wrote " << total_written << " filtered records to: "
              << output_filename << std::endl;
    return true;
}

// ============================================================================
// BamRouter Implementation
// ============================================================================

BamRouter::BamRouter()
    : remove_duplicates_(false) {
}

BamRouter::~BamRouter() {
}

bool BamRouter::routeBamSinglePass(const std::string& input_bam,
                                    const std::map<CellTypeDonorKey, BamGroupInfo>& groups,
                                    const std::string& output_dir) {
    
    std::cout << "\n=== SCAR: Processing BAM file: " << input_bam << " ===" << std::endl;
    
    samFile* in_bam = sam_open(input_bam.c_str(), "r");
    if (!in_bam) {
        std::cerr << "Error: Cannot open input BAM: " << input_bam << std::endl;
        return false;
    }
    
    sam_hdr_t* header = sam_hdr_read(in_bam);
    if (!header) {
        std::cerr << "Error: Cannot read BAM header" << std::endl;
        sam_close(in_bam);
        return false;
    }
    
    fs::create_directories(output_dir);
    
    std::unordered_map<std::string, CellTypeDonorKey> barcode_to_key;
    for (const auto& [key, group] : groups) {
        if (group.bamfile != input_bam) continue;
        for (const auto& barcode : group.barcodes) {
            barcode_to_key[barcode] = key;
        }
    }
    
    std::cout << "Loaded " << barcode_to_key.size() << " barcodes for this BAM file" << std::endl;
    
    std::map<std::string, std::shared_ptr<BedRegions>> celltype_beds;
    std::set<std::string> unique_celltypes;
    for (const auto& [key, group] : groups) {
        unique_celltypes.insert(group.celltype);
    }
    
    for (const auto& celltype : unique_celltypes) {
        std::string bed_file = "";
        for (const auto& [key, group] : groups) {
            if (group.celltype == celltype && !group.bed_file.empty()) {
                bed_file = group.bed_file;
                break;
            }
        }
        
        if (!bed_file.empty()) {
            auto bed_regions = std::make_shared<BedRegions>();
            if (bed_regions->loadBed(bed_file)) {
                celltype_beds[celltype] = bed_regions;
                std::cout << "Loaded BED for " << celltype << ": " 
                          << bed_regions->getIntervalCount() << " intervals" << std::endl;
            }
        }
    }
    
    std::map<CellTypeDonorKey, OutputFile> output_files;
    for (const auto& [key, group] : groups) {
        if (group.bamfile != input_bam) continue;
        
        std::string celltype_dir = output_dir + "/" + group.celltype;
        fs::create_directories(celltype_dir);
        
        std::string output_filename = celltype_dir + "/" + group.donor_id + ".bam";
        samFile* out_bam = sam_open(output_filename.c_str(), "wb");
        if (!out_bam) {
            std::cerr << "Error: Cannot create output BAM: " << output_filename << std::endl;
            for (auto& [k, of] : output_files) {
                sam_close(of.handle);
            }
            sam_hdr_destroy(header);
            sam_close(in_bam);
            return false;
        }
        
        if (sam_hdr_write(out_bam, header) < 0) {
            std::cerr << "Error: Cannot write BAM header to " << output_filename << std::endl;
            sam_close(out_bam);
            for (auto& [k, of] : output_files) {
                sam_close(of.handle);
            }
            sam_hdr_destroy(header);
            sam_close(in_bam);
            return false;
        }
        
        auto bed_it = celltype_beds.find(group.celltype);
        std::shared_ptr<BedRegions> bed_regions = (bed_it != celltype_beds.end()) ? bed_it->second : nullptr;
        
        output_files[key] = {output_filename, out_bam, {}, bed_regions, 0};
        output_files[key].write_buffer.reserve(OutputFile::BUFFER_SIZE);
        std::cout << "Opened output: " << output_filename 
                  << " (expecting " << group.cell_count << " cells)"
                  << (bed_regions ? " [with BED filtering]" : "") << std::endl;
    }
    
    if (remove_duplicates_) {
        std::cout << "Filtering: Removing duplicate reads (flag 0x400)" << std::endl;
    }
    
    std::cout << "\n=== SCAR: Routing reads with per-celltype peak filtering ===" << std::endl << std::endl;
    
    bam1_t* aln = bam_init1();
    size_t total_reads = 0;
    size_t filtered_duplicate = 0;
    size_t filtered_bed = 0;
    size_t no_barcode = 0;
    size_t barcode_not_in_list = 0;
    size_t written = 0;
    size_t buffer_flushes = 0;
    
    auto start_time = std::chrono::high_resolution_clock::now();
    auto last_checkpoint_time = start_time;
    const int CHECKPOINT_INTERVAL_MS = 5000;
    
    while (sam_read1(in_bam, header, aln) >= 0) {
        total_reads++;
        
        uint8_t* cb_tag = bam_aux_get(aln, "CB");
        if (!cb_tag) {
            no_barcode++;
            continue;
        }
        
        std::string cb_value = bam_aux2Z(cb_tag);
        
        auto it = barcode_to_key.find(cb_value);
        if (it == barcode_to_key.end()) {
            barcode_not_in_list++;
            continue;
        }
        
        if (remove_duplicates_ && (aln->core.flag & BAM_FDUP)) {
            filtered_duplicate++;
            continue;
        }
        
        const CellTypeDonorKey& key = it->second;
        auto& output = output_files[key];
        
        if (output.bed_regions) {
            int32_t tid = aln->core.tid;
            if (tid >= 0 && tid < header->n_targets) {
                if (!output.bed_regions->overlaps(sam_hdr_tid2name(header, tid), aln->core.pos, aln->core.l_qseq)) {
                    filtered_bed++;
                    continue;
                }
            } else {
                filtered_bed++;
                continue;
            }
        }
        
        output.write_buffer.push_back(bam_dup1(aln));
        written++;
        
        if (output.write_buffer.size() >= OutputFile::BUFFER_SIZE) {
            if (!output.flush(header)) {
                std::cerr << "Error: Failed to write alignments to " << output.filename << std::endl;
                bam_destroy1(aln);
                for (auto& [k, of] : output_files) {
                    for (auto* b : of.write_buffer) {
                        bam_destroy1(b);
                    }
                    sam_close(of.handle);
                }
                sam_hdr_destroy(header);
                sam_close(in_bam);
                return false;
            }
            buffer_flushes++;
        }
        
        auto now = std::chrono::high_resolution_clock::now();
        auto elapsed_since_checkpoint = std::chrono::duration_cast<std::chrono::milliseconds>(now - last_checkpoint_time);
        
        if (elapsed_since_checkpoint.count() >= CHECKPOINT_INTERVAL_MS) {
            auto total_elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - start_time);
            double reads_per_second = total_elapsed.count() > 0 ? 
                total_reads / static_cast<double>(total_elapsed.count()) : 0;
            
            std::cout << "Processed " << total_reads << " reads | Written " << written
                      << " | Dup " << filtered_duplicate 
                      << " | BED " << filtered_bed
                      << " | No CB " << no_barcode
                      << " | Not in list " << barcode_not_in_list
                      << "\n  Throughput: " << std::fixed << std::setprecision(2)
                      << reads_per_second << " reads/sec (avg)"
                      << " | Flushes: " << buffer_flushes
                      << " | Elapsed: " << total_elapsed.count() << "s" << std::endl;
            
            last_checkpoint_time = now;
        }
    }
    
    bam_destroy1(aln);
    sam_close(in_bam);
    
    std::cout << "\n=== Flushing remaining buffers ===" << std::endl;
    size_t total_buffered = 0;
    for (auto& [key, output] : output_files) {
        total_buffered += output.write_buffer.size();
        if (!output.flush(header)) {
            std::cerr << "Error: Failed to flush final buffer for " << output.filename << std::endl;
            sam_hdr_destroy(header);
            return false;
        }
    }
    std::cout << "Flushed " << total_buffered << " buffered reads" << std::endl;
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto total_elapsed = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
    double throughput = total_elapsed.count() > 0 ? total_reads / static_cast<double>(total_elapsed.count()) : 0;
    
    std::cout << "\n=== Processing Summary ===" << std::endl;
    std::cout << "Total reads processed: " << total_reads << std::endl;
    std::cout << "  Written: " << written << std::endl;
    std::cout << "  Filtered (duplicates): " << filtered_duplicate << std::endl;
    std::cout << "  Filtered (BED): " << filtered_bed << std::endl;
    std::cout << "  No CB tag: " << no_barcode << std::endl;
    std::cout << "  CB not in list: " << barcode_not_in_list << std::endl;
    std::cout << "  Buffer flushes: " << buffer_flushes << std::endl;
    std::cout << "Total time: " << total_elapsed.count() << "s" << std::endl;
    std::cout << "Overall throughput: " << std::fixed << std::setprecision(2) 
              << throughput << " reads/sec" << std::endl;
    
    std::cout << "\n=== Closing output files and creating indexes ===" << std::endl;
    for (auto& [key, output] : output_files) {
        std::cout << "Closing: " << key.celltype << "/" << key.donor_id
                  << ".bam (" << output.written_count << " reads written)" << std::endl;
        sam_close(output.handle);
        
        std::cout << "  Creating index..." << std::flush;
        if (sam_index_build(output.filename.c_str(), 0) < 0) {
            std::cerr << " Warning: Failed to create index for " << output.filename << std::endl;
        } else {
            std::cout << " done" << std::endl;
        }
    }
    
    sam_hdr_destroy(header);
    
    std::cout << "\n=== SCAR: BAM file processed successfully ===" << std::endl;
    return true;
}

