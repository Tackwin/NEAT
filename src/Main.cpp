
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

void start_dockspace() noexcept;
void end_dockspace() noexcept;

int main(int, char**) {
	// Setup window
	glfwSetErrorCallback(glfw_error_callback);
	if (!glfwInit()) return 1;

	GLFWwindow* window = glfwCreateWindow(1280, 720, "NEAT", NULL, NULL);
	if (window == NULL) return 1;

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
			for (size_t i = neat.population.size(); i < neat.population_size; ++i) {
				neat.population.push_back(neat_win.initial_genome);
			}

			neat.evaluate([](Network& n) {
				float inputs[] = {
					1, 0, 0,
					1, 1, 0,
					1, 0, 1,
					1, 1, 1
				};
				float target_outputs[] = {
					0, 1, 1, 0
				};
				float predicted_outputs[4];

				n.compute(inputs + 0, 3, predicted_outputs + 0, 1);
				n.compute(inputs + 3, 3, predicted_outputs + 1, 1);
				n.compute(inputs + 6, 3, predicted_outputs + 2, 1);
				n.compute(inputs + 9, 3, predicted_outputs + 3, 1);

				float s = 0;
				for (size_t i = 0; i < 4; ++i) {
					float dt = (target_outputs[i] - predicted_outputs[i]);
					s += dt * dt;
				}

				return 1.f / s;
			});

			neat.select();
			neat.populate();

			float m = 0;
			for (auto& x : neat.results) if (m < x.fitness) m = x.fitness;
			neat_win.max_fitness.push_back(m);
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
