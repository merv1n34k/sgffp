#pragma description SnapGene binary format
#pragma endian big
#pragma MIME application/x-snapgene

#include <std/sys.pat>
#include <std/mem.pat>
#include <std/io.pat>

// Color coding for different block types
#define COLOR_HEADER 0xFF8888
#define COLOR_DATA 0x88FF88
#define COLOR_META 0x8888FF
#define COLOR_UNKNOWN 0xFFFF88

// Block type enumeration
enum BlockType : u8 {
    SEQ_DNA = 0x00,
    SEQ_RNA = 0x20,
    SEQ_PROT = 0x15,
    COMPRESSED_SEQ = 0x01,
    UNKNOWN_2 = 0x02,
    UNKNOWN_3 = 0x03,
    PRIMERS = 0x05,
    NOTES = 0x06,
    HISTORY_TREE = 0x07,
    ADDITIONAL_PROPS = 0x08,
    FILE_DESCRIPTION = 0x09,
    FEATURES = 0x0A,
    HISTORY_NODE = 0x0B,
    UNKNOWN_13 = 0x0D,
    ALIGN_SEQ_16 = 0x10,
    ALIGN_SEQ_17 = 0x11,
    SEQ_TRACE = 0x12,
    URACIL_POSITIONS = 0x13,
    CUSTOM_COLORS = 0x14,
    ENZ_VIZ = 0x1c,
    ENZ_CUST = 0xe
};

// Header structure
struct SnapGeneHeader {
    u8 magic [[color("FF0000")]];  // Should be 0x09 ('\t')

    // Validate magic byte
    if (magic != 0x09) {
        std::error("Invalid magic byte! Expected 0x09, got 0x{:02X}", magic);
    }

    u32 headerLength [[color("00FF00")]];  // Should be 14
    char title[8] [[color("0000FF")]];     // Should be "TestView"

    // Validate header
    if (headerLength != 14) {
        std::warning("Unexpected header length: {}", headerLength);
    }
    if (title != "SnapGene") {
        std::error("Invalid title! Expected 'SnapGene'");
    }

    u16 sequenceType [[color("FFFF00")]];
    u16 exportVersion [[color("FF00FF")]];
    u16 importVersion [[color("00FFFF")]];
} [[static, color(COLOR_HEADER)]];

// Generic block structure
struct Block {
    BlockType type [[color("FF8800")]];
    u32 dataSize [[color("0088FF")]];

    // Color based on block type
    if (type == BlockType::SEQ_DNA) {
        u8 data[dataSize] [[color("FFFFFF")]];
    } else if (type == BlockType::SEQ_RNA) {
        u8 data[dataSize] [[color("888888")]];
    } else if (type == BlockType::SEQ_PROT) {
        u8 data[dataSize] [[color("FF8888")]];
    } else if (type == BlockType::COMPRESSED_SEQ) {
        u8 data[dataSize] [[color("FF0000")]];
    } else if (type == BlockType::PRIMERS) {
        char stringData[dataSize] [[color("88FF88")]];
    } else if (type == BlockType::NOTES) {
        char notes[dataSize] [[color("88FFFF")]];
    } else if (type == BlockType::ALIGN_SEQ_16 || type == BlockType::ALIGN_SEQ_17) {
        u8 rawData[dataSize] [[color("8888FF")]];
    } else if (type == BlockType::FEATURES) {
        u8 featureData[dataSize] [[color("FF88FF")]];
    } else if (type == BlockType::FILE_DESCRIPTION) {
        char description[dataSize] [[color("FFFF88")]];
    } else if (type == BlockType::SEQ_TRACE) {
        u8 traceData[dataSize] [[color("88FF88")]];
    } else if (type == BlockType::URACIL_POSITIONS) {
        u8 pointData[dataSize] [[color("FF8888")]];
    } else {
        u8 unknownData[dataSize] [[color(COLOR_UNKNOWN)]];
    }
} [[format("format_block")]];

// Format function for nice display
fn format_block(Block block) {
    return std::format("Type: {} (0x{:02X}), Size: {} bytes",
        get_block_name(block.type), u8(block.type), block.dataSize);
};

// Get human-readable block name
fn get_block_name(BlockType type) {
    if (type == BlockType::SEQ_DNA) return "DNA Sequence";
    if (type == BlockType::SEQ_RNA) return "RNA Sequence";
    if (type == BlockType::SEQ_PROT) return "Protein Sequence";
    if (type == BlockType::COMPRESSED_SEQ) return "Compressed sequence";
    if (type == BlockType::PRIMERS) return "Primers";
    if (type == BlockType::NOTES) return "Notes";
    if (type == BlockType::HISTORY_TREE) return "History Tree";
    if (type == BlockType::ADDITIONAL_PROPS) return "Additional Properties";
    if (type == BlockType::FILE_DESCRIPTION) return "File Description";
    if (type == BlockType::FEATURES) return "Features";
    if (type == BlockType::HISTORY_NODE) return "History Node";
    if (type == BlockType::ALIGN_SEQ_16) return "Alignable sequence";
    if (type == BlockType::ALIGN_SEQ_17) return "Alignable sequence";
    if (type == BlockType::SEQ_TRACE) return "Sequence Trace";
    if (type == BlockType::URACIL_POSITIONS) return "Uracil Positions";
    if (type == BlockType::CUSTOM_COLORS) return "Custom Colors";
    if (type == BlockType::ENZ_VIZ) return "Enzyme vizualisation";
    if (type == BlockType::ENZ_CUST) return "Custom enzymes";
    return std::format("Unknown (0x{:02X})", u8(type));
};

// Main file structure
struct SnapGeneFile {
    SnapGeneHeader header;
    Block blocks[while(!std::mem::eof())];
} [[inline]];

// Create the main structure
SnapGeneFile file @ 0x00;
