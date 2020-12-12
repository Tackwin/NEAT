
#include "AI/Neat.hpp"
#include "UI/Window.hpp"

#include "imgui/imgui.h"
#include "imgui/imnodes.h"
#include "imgui/implot.h"
#include "imgui/imgui_impl_glfw.h"
#include "imgui/imgui_impl_opengl2.h"

#include "Argh.hpp"
#include "xstd.hpp"

#include "Problem/Demo.hpp"

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


void gui(argh::parser& flags) noexcept;
void cli(argh::parser& flags) noexcept;

int main(int, char** argv) {
	argh::parser flags(argv + 1);

	if (flags["no-ui"]) cli(flags);
	else                gui(flags);

	return 0;
}

void cli(argh::parser& flags) noexcept {
	if (!flags("init-gen")) {
		printf("You need to provide an input network");
		return;
	}

	Neat neat;
	size_t n_iter = 0;
	int eval = 0;
	float goal = 0.f;
	std::string init_genome_path;
	
	flags("init-gen") >> init_genome_path;
	flags("pop") >> neat.population_size;
	flags("goal", 0) >> goal;
	flags("n", 0) >> n_iter;
	flags("eval", 0) >> eval;

	Neat_Window::Evaluation evaluation = (Neat_Window::Evaluation)eval;
	auto fit = xor_fitness;
	switch (evaluation) {
	case Neat_Window::Evaluation::XOR :
		fit = xor_fitness;
		break;
	case Neat_Window::Evaluation::SPV :
		fit = spv_fitness;
		break;
	default:
		fit = xor_fitness;
		break;
	}

	FILE* f = fopen(init_genome_path.c_str(), "rb");

	fseek(f, 0, SEEK_END);
	long fsize = ftell(f);
	fseek(f, 0, SEEK_SET);  /* same as rewind(f); */

	std::vector<std::uint8_t> data;
	data.resize(fsize);
	fread(data.data(), 1, fsize, f);
	fclose(f);

	Genome init_genome = Genome::deserialize(data);

	auto start = xstd::seconds();
	neat.complete_with(init_genome, neat.population_size);
	for (size_t i = 0; n_iter == 0 || i < n_iter; ++i) {
		if (n_iter <= 10 || (i % (n_iter / 10)) == 0) {
			printf(
				"Doing Generation % *d, with % 5d species\n",
				(int)(std::ceil(std::log10(n_iter + 1)) + 1),
				(int)(neat.generation_number + 1),
				(int)(neat.species.size())
			);
		}

		neat.evaluate(fit, nullptr);
		neat.select();
		neat.populate();
		neat.speciate();
		neat.generation_number++;

		auto m = std::max_element(
			std::begin(neat.results),
			std::end(neat.results),
			[](auto a, auto b) { return a.fitness < b.fitness; }
		);
		if (goal > 0 && goal < m->fitness) {
			printf("After %zu generations, reached %lf\n", neat.generation_number, m->fitness);
			printf("Saved to best_genome\n");

			data = m->g->serialize();
			FILE* f = fopen("best_genome", "wb");
			fwrite(data.data(), 1, data.size(), f);
			fclose(f);

			break;
		}
	}
	auto end = xstd::seconds();

	auto t = end - start;
	auto kgs = (neat.population_size * neat.generation_number / t) / 1'000.0;

	printf(
		"Done %zu generation of %zu genomesin:\n",
		neat.generation_number,
		neat.population_size
	);
	printf("%10.5f seconds, %10.5f Kgs (Kilo genomes per seconds)\n", t, kgs);
}

void gui(argh::parser& flags) noexcept {
	// Setup window
	glfwSetErrorCallback(glfw_error_callback);
	if (!glfwInit()) return;

	GLFWwindow* window = glfwCreateWindow(1600, 900, "NEAT", NULL, NULL);
	if (window == NULL) return;

	glfwSetKeyCallback(window, key_callback);
	glfwMakeContextCurrent(window);
	glfwSwapInterval(1); // Enable vsync

	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	ImPlot::CreateContext();
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

		auto f = xor_fitness;
		switch (neat_win.evaluation) {
		case Neat_Window::Evaluation::XOR :
			f = xor_fitness;
			break;
		case Neat_Window::Evaluation::SPV :
			f = spv_fitness;
			break;
		default:
			f = xor_fitness;
			break;
		}

		switch (neat_win.evaluation) {
		case Neat_Window::Evaluation::SPV :
			neat_win.render_phenotype = spv_render;
			break;
		default:
			neat_win.render_phenotype = nullptr;
			break;
		}

		if (neat_win.run_generation || neat_win.auto_run) {
			neat.complete_with(neat_win.initial_genome, neat.population_size);

			neat.evaluate(f, nullptr);
			neat.speciate();
			neat_win.get_stats(neat, neat.results);
			neat.select();
			neat.populate();


			neat.generation_number++;
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
	ImPlot::DestroyContext();
	ImGui::DestroyContext();

	glfwDestroyWindow(window);
	glfwTerminate();
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
