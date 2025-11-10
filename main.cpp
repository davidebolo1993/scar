#include "scar.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <map>
#include <set>

void print_usage(const char* prog_name) {
    std::cout << "SCAR - Single-Cell ATAC Router\n";
    std::cout << "Version 1.0\n\n";
    std::cout << "Usage: " << prog_name << " [options]\n\n";
    std::cout << "Subcommands:\n";
    std::cout << "  filter  Filter metadata and prepare for splitting\n";
    std::cout << "  split   Split BAM files by celltype and donor\n\n";
    
    std::cout << "Filter options:\n";
    std::cout << "  -i, --input FILE           Input metadata TSV file\n";
    std::cout << "  -o, --output FILE          Output filtered metadata file\n";
    std::cout << "  -c, --celltype COLUMN      Celltype column (celltypes_LEV1 or celltypes_LEV2) [default: celltypes_LEV1]\n";
    std::cout << "  -m, --min-cells INT        Minimum cells per donor/celltype [default: 10]\n";
    std::cout << "  -p, --peak-matrix FILE     Peak matrix TSV file (celltype columns)\n";
    std::cout << "  -b, --bed-output-dir DIR   Output directory for celltype-specific BED files\n\n";
    
    std::cout << "Split options:\n";
    std::cout << "  -i, --input FILE              Filtered metadata file from filter step\n";
    std::cout << "  -o, --output-dir DIR          Output directory for split BAM files\n";
    std::cout << "  -c, --celltype COLUMN         Celltype column (celltypes_LEV1 or celltypes_LEV2) [default: celltypes_LEV1]\n";
    std::cout << "  -d, --remove-duplicates       Remove duplicate reads (samtools -F 1024 equivalent)\n";
}

int run_filter(int argc, char* argv[]) {
    std::string input_file;
    std::string output_file;
    std::string celltype_column = "celltypes_LEV1";
    std::string peak_matrix_file;
    std::string bed_output_dir;
    int min_cells = 10;
    
    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "-i") == 0 || strcmp(argv[i], "--input") == 0) {
            input_file = argv[++i];
        } else if (strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--output") == 0) {
            output_file = argv[++i];
        } else if (strcmp(argv[i], "-c") == 0 || strcmp(argv[i], "--celltype") == 0) {
            celltype_column = argv[++i];
        } else if (strcmp(argv[i], "-m") == 0 || strcmp(argv[i], "--min-cells") == 0) {
            min_cells = std::stoi(argv[++i]);
        } else if (strcmp(argv[i], "-p") == 0 || strcmp(argv[i], "--peak-matrix") == 0) {
            peak_matrix_file = argv[++i];
        } else if (strcmp(argv[i], "-b") == 0 || strcmp(argv[i], "--bed-output-dir") == 0) {
            bed_output_dir = argv[++i];
        }
    }
    
    if (input_file.empty() || output_file.empty()) {
        std::cerr << "Error: Input and output files are required" << std::endl;
        return 1;
    }
    
    MetadataProcessor processor(min_cells);
    
    if (!processor.loadMetadata(input_file, celltype_column)) {
        return 1;
    }
    
    if (!peak_matrix_file.empty() && !bed_output_dir.empty()) {
        if (!processor.loadPeakMatrix(peak_matrix_file, bed_output_dir)) {
            return 1;
        }
    }
    
    processor.filterMetadata();
    
    if (!processor.writeFilteredMetadata(output_file)) {
        return 1;
    }
    
    return 0;
}

int run_split(int argc, char* argv[]) {
    std::string input_file;
    std::string output_dir;
    std::string celltype_column = "celltypes_LEV1";
    bool remove_duplicates = false;
    
    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "-i") == 0 || strcmp(argv[i], "--input") == 0) {
            input_file = argv[++i];
        } else if (strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--output-dir") == 0) {
            output_dir = argv[++i];
        } else if (strcmp(argv[i], "-c") == 0 || strcmp(argv[i], "--celltype") == 0) {
            celltype_column = argv[++i];
        } else if (strcmp(argv[i], "-d") == 0 || strcmp(argv[i], "--remove-duplicates") == 0) {
            remove_duplicates = true;
        }
    }
    
    if (input_file.empty() || output_dir.empty()) {
        std::cerr << "Error: Input file and output directory are required" << std::endl;
        return 1;
    }
    
    std::ifstream infile(input_file);
    if (!infile.is_open()) {
        std::cerr << "Error: Cannot open filtered metadata file" << std::endl;
        return 1;
    }
    
    std::string line;
    bool is_header = true;
    std::map<CellTypeDonorKey, BamGroupInfo> groups;
    std::set<std::string> unique_bams;
    
    while (std::getline(infile, line)) {
        if (is_header) {
            is_header = false;
            continue;
        }
        
        std::istringstream iss(line);
        std::string barcode, pool_id, donor_id, celltype, bamfile, bed_file;
        
        if (!(iss >> barcode >> pool_id >> donor_id >> celltype >> bamfile >> bed_file)) {
            continue;
        }
        
        CellTypeDonorKey key{celltype, donor_id};
        
        if (groups.find(key) == groups.end()) {
            groups[key] = BamGroupInfo{bamfile, donor_id, celltype, bed_file, {}, 0};
        }
        
        groups[key].barcodes.insert(barcode);
        groups[key].cell_count++;
        unique_bams.insert(bamfile);
    }
    
    std::cout << "SCAR: Loaded " << groups.size() << " donor/celltype combinations\n";
    std::cout << "SCAR: Processing " << unique_bams.size() << " unique BAM files\n\n";
    
    BamRouter router;
    router.setRemoveDuplicates(remove_duplicates);
    
    for (const auto& bam_file : unique_bams) {
        if (!router.routeBamSinglePass(bam_file, groups, output_dir)) {
            std::cerr << "Error processing BAM: " << bam_file << std::endl;
            return 1;
        }
    }
    
    std::cout << "\n=== SCAR: All BAM files processed successfully ===" << std::endl;
    return 0;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }
    
    std::string subcommand = argv[1];
    
    if (subcommand == "filter") {
        return run_filter(argc, argv);
    } else if (subcommand == "split") {
        return run_split(argc, argv);
    } else {
        std::cerr << "Error: Unknown subcommand '" << subcommand << "'" << std::endl;
        print_usage(argv[0]);
        return 1;
    }
}

