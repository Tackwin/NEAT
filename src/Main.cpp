
#include "AI/Neat.hpp"
#include "UI/Window.hpp"

#include "imgui/imgui.h"
#include "imgui/imnodes.h"
#include "imgui/imgui_impl_glfw.h"
#include "imgui/imgui_impl_opengl2.h"

#include <stdio.h>

#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#endif
#include "GLFW/glfw3.h"

static void glfw_error_callback(int error, const char* description) {
	fprintf(stderr, "Glfw Error %d: %s\n", error, description);
}

GLFWmonitor* get_current_monitor(GLFWwindow *window) noexcept;

void toggle_fullscreen(GLFWwindow* w) noexcept {
	thread_local bool fullscreen = false;
	thread_local int last_x = 0;
	thread_local int last_y = 0;
	thread_local int last_w = 0;
	thread_local int last_h = 0;

	if (!fullscreen) {
		fullscreen = true;
		glfwGetWindowPos(w, &last_x, &last_y);
		glfwGetWindowSize(w, &last_w, &last_h);

		auto m = get_current_monitor(w);
		auto mode = glfwGetVideoMode(m);
		glfwSetWindowMonitor(
			w,
			m,
			0,
			0,
			mode->width,
			mode->height,
			GLFW_DONT_CARE
		);
	} else {
		fullscreen = false;
		glfwSetWindowMonitor(
			w,
			nullptr,
			last_x,
			last_y,
			last_w,
			last_h,
			GLFW_DONT_CARE
		);

	}

}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
	if (key == GLFW_KEY_F11 && action == GLFW_PRESS) toggle_fullscreen(window);
}

int main(int, char**) {
	// Setup window

	glfwSetErrorCallback(glfw_error_callback);
	if (!glfwInit()) return 1;

	GLFWwindow* window = glfwCreateWindow(1600, 900, "NEAT", NULL, NULL);
	if (window == NULL) return 1;

	glfwSetKeyCallback(window, key_callback);
	glfwMakeContextCurrent(window);
	glfwSwapInterval(1); // Enable vsync

	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	ImGuiIO& io = ImGui::GetIO();
	io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;
	
	imnodes::Initialize();

	Neat neat;
	Neat_Window neat_win;

	ImGui::StyleColorsDark();

	// Setup Platform/Renderer backends
	ImGui_ImplGlfw_InitForOpenGL(window, true);
	ImGui_ImplOpenGL2_Init();

	size_t i = 0;

	// Our state
	bool show_demo_window = false;
	ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);

	// Main loop
	while (!glfwWindowShouldClose(window)) {
		glfwPollEvents();

		// Start the Dear ImGui frame
		ImGui_ImplOpenGL2_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();
		
		ImGui::DockSpaceOverViewport(
			ImGui::GetMainViewport(),
			ImGuiDockNodeFlags_PassthruCentralNode
		);
		if (show_demo_window) ImGui::ShowDemoWindow(&show_demo_window);

		neat_win.render(neat);
		if (neat_win.run_generation || neat_win.auto_run) {
			neat.complete_with(neat_win.initial_genome, neat.population_size);
			neat.evaluate(xor_fitness);
			neat.select();
			neat.populate();
			neat.speciate();
			neat.generation_number++;
			neat_win.get_stats(neat.results);
		}

		// Rendering
		ImGui::Render();
		int display_w, display_h;
		glfwGetFramebufferSize(window, &display_w, &display_h);
		glViewport(0, 0, display_w, display_h);
		glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);
		glClear(GL_COLOR_BUFFER_BIT);

		ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());

		glfwMakeContextCurrent(window);
		glfwSwapBuffers(window);
	}

	// Cleanup
	ImGui_ImplOpenGL2_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();

	glfwDestroyWindow(window);
	glfwTerminate();

	return 0;
}

GLFWmonitor* get_current_monitor(GLFWwindow *window) noexcept {
	int nmonitors, i;
	int wx, wy, ww, wh;
	int mx, my, mw, mh;
	int overlap, bestoverlap;
	GLFWmonitor *bestmonitor;
	GLFWmonitor **monitors;
	const GLFWvidmode *mode;

	bestoverlap = 0;
	bestmonitor = NULL;

	glfwGetWindowPos(window, &wx, &wy);
	glfwGetWindowSize(window, &ww, &wh);
	monitors = glfwGetMonitors(&nmonitors);

	for (i = 0; i < nmonitors; i++) {
		mode = glfwGetVideoMode(monitors[i]);
		glfwGetMonitorPos(monitors[i], &mx, &my);
		mw = mode->width;
		mh = mode->height;

		overlap =
			std::max(0, std::min(wx + ww, mx + mw) - std::max(wx, mx)) *
			std::max(0, std::min(wy + wh, my + mh) - std::max(wy, my));

		if (bestoverlap < overlap) {
			bestoverlap = overlap;
			bestmonitor = monitors[i];
		}
	}

	return bestmonitor;
}
