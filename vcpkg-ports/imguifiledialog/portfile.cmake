vcpkg_check_linkage(ONLY_STATIC_LIBRARY)

vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO aiekick/ImGuiFileDialog
    REF master
    SHA512 039fd70b784633d60b7c8ab2893efb540195bcb1443d98dba47d6350d7e96a77b80f45a40cf046465a34b7a6fd914c7d47af7f7080257495f4b826d2349729be
    HEAD_REF master
)

file(WRITE "${SOURCE_PATH}/ImGuiFileDialogConfig.h" "#pragma once\n#define USE_STD_FILESYSTEM\n")

vcpkg_cmake_configure(
    SOURCE_PATH "${SOURCE_PATH}"
    OPTIONS
        -DBUILD_SHARED_LIBS=OFF
        -DIGFD_INSTALL=ON
        -DCMAKE_CXX_STANDARD=17
)

vcpkg_cmake_install()
vcpkg_copy_pdbs()

vcpkg_cmake_config_fixup(CONFIG_PATH lib/cmake/ImGuiFileDialog PACKAGE_NAME ImGuiFileDialog)
file(READ "${CURRENT_PACKAGES_DIR}/share/imguifiledialog/ImGuiFileDialogConfig.cmake" cmake_config)
string(PREPEND cmake_config [[
include(CMakeFindDependencyMacro)
find_dependency(imgui CONFIG REQUIRED)
]])
file(WRITE "${CURRENT_PACKAGES_DIR}/share/imguifiledialog/ImGuiFileDialogConfig.cmake" "${cmake_config}")

file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/include")
file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/bin" "${CURRENT_PACKAGES_DIR}/debug/bin")

vcpkg_install_copyright(FILE_LIST "${SOURCE_PATH}/LICENSE")
