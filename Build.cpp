#include "Ease.hpp"

Build build(Flags flags) noexcept {
	auto b = Build::get_default(flags);

	b.name = "NEAT";

	b.add_define("_CRT_SECURE_NO_WARNINGS");
	b.add_define("IMGUI_DISABLE_INCLUDE_IMCONFIG_H");

	b.add_header("src/");
	b.add_header("src/imgui");

	b.add_source_recursively("src/");

	if (Env::Win32) {
		b.add_library("opengl32");
		b.add_library("gdi32");
		b.add_library("advapi32");
		b.add_library("winmm");
		b.add_library("msvcrt");
		b.add_library("shell32");
		b.add_library("lib/win32/glfw3");
	} else {
		b.add_library("glfw");
		b.add_library("GL");
	}

	return b;
}