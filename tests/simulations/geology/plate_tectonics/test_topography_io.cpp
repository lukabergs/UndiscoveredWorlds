#include "gtest/gtest.h"
#include "heightmap_io.hpp"

#include <filesystem>
#include <vector>

namespace fs = std::filesystem;

namespace {

struct ScopedDelete {
    explicit ScopedDelete(fs::path _path) : path(std::move(_path)) {}
    ~ScopedDelete()
    {
        std::error_code ec;
        fs::remove(path, ec);
    }

    fs::path path;
};

} // namespace

TEST(TopographyIo, Tiff16RoundTrip)
{
    const int width = 3;
    const int height = 2;
    const std::vector<uint16_t> samples = {0U, 1U, 25U, 1024U, 32768U, 65535U};
    const fs::path path = fs::current_path() / "topography_roundtrip.tiff";
    ScopedDelete cleanup(path);

    ASSERT_EQ(0, writeImageGrayTiff16(path.string().c_str(), width, height, samples.data()));

    std::vector<uint16_t> loaded;
    int loaded_width = 0;
    int loaded_height = 0;
    ASSERT_EQ(0, readImageGrayTiff16(path.string().c_str(), loaded, loaded_width, loaded_height));
    EXPECT_EQ(width, loaded_width);
    EXPECT_EQ(height, loaded_height);
    EXPECT_EQ(samples, loaded);
}

TEST(TopographyIo, MetadataRoundTrip)
{
    const fs::path path = fs::current_path() / "topography_roundtrip.tiff.json";
    ScopedDelete cleanup(path);

    TopographyCodec::Metadata metadata;
    metadata.width = 64;
    metadata.height = 32;
    metadata.sea_level_m = 31000;
    metadata.max_height_m = TopographyCodec::kMaxHeightMeters;
    metadata.version = TopographyCodec::kMetadataVersion;
    metadata.format = TopographyCodec::kMetricFormatTiff16;
    metadata.endianness = TopographyCodec::kLittleEndian;
    metadata.layout = TopographyCodec::kRowMajorLayout;

    ASSERT_EQ(0, writeTopographyMetadataJson(path.string().c_str(), metadata));

    TopographyCodec::Metadata loaded;
    ASSERT_EQ(0, readTopographyMetadataJson(path.string().c_str(), loaded));
    EXPECT_EQ(metadata.width, loaded.width);
    EXPECT_EQ(metadata.height, loaded.height);
    EXPECT_EQ(metadata.sea_level_m, loaded.sea_level_m);
    EXPECT_EQ(metadata.max_height_m, loaded.max_height_m);
    EXPECT_EQ(metadata.version, loaded.version);
    EXPECT_EQ(metadata.format, loaded.format);
    EXPECT_EQ(metadata.endianness, loaded.endianness);
    EXPECT_EQ(metadata.layout, loaded.layout);
}
