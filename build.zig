const std = @import("std");

pub fn build(b: *std.Build) void {
    _ = b;
    // const target = b.standardTargetOptions(.{});
    // const optimize = b.standardOptimizeOption(.{});
    // const cpp_flags = &.{
    //     "-std=c++11",
    //     "-fno-rtti",
    // };

    // const module = b.addModule("signalsmith_stretch", .{
    //     .root_source_file = b.path("zig-src/root.zig"),
    //     .target = target,
    //     .optimize = optimize,
    //     .link_libc = true,
    //     .link_libcpp = true,
    // });

    // const audioeffectx_h = b.path("zig-src/shim/");
    // const audioeffectx_cpp = b.path("zig-src/shim/audioeffect.cpp");
    // module.addIncludePath(audioeffectx_h.dirname());
    // module.addCSourceFile(.{
    //     .file = audioeffectx_cpp,
    //     .language = .cpp,
    //     .flags = cpp_flags,
    // });

}
