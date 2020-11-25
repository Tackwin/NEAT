#include "Window.hpp"

#include "imgui.h"
#include "imgui/imnodes.h"

void Network_Window::embed_render(Network& network) noexcept {
	links_created.clear();
	links_deleted.clear();
	links_hovered.clear();

	imnodes::BeginNodeEditor();

	for (size_t i = 0; i < network.nodes.size(); ++i) {
		imnodes::BeginNode(i);

		imnodes::BeginNodeTitleBar();
		ImGui::Text(
			"%s",
			network.node_kinds[i] == Network::Node_Kind::Hidden ? "Hidden" :
			(network.node_kinds[i] == Network::Node_Kind::Output ? "Output" :
			"Input")
		);
		imnodes::EndNodeTitleBar();

		ImGui::Text("%zu", i);

		if (network.node_kinds[i] != Network::Node_Kind::Output) {
			imnodes::BeginOutputAttribute(i * 2 + 0);
			imnodes::EndOutputAttribute();
		}
		if (network.node_kinds[i] != Network::Node_Kind::Sensor) {
			imnodes::BeginInputAttribute(i * 2 + 1);
			imnodes::EndInputAttribute();
		}

		imnodes::EndNode();
	}

	for (size_t i = 0; i < network.connections.size(); ++i) {
		imnodes::Link(
			i,
			network.connections[i].in  * 2 + 0,
			network.connections[i].out * 2 + 1
		);
	}

	imnodes::EndNodeEditor();

	int start_id = 0;
	int end_id = 0;
	int hover_id = 0;
	if (imnodes::IsLinkCreated(&start_id, &end_id)) {
		links_created.push_back(start_id / 2);
		links_created.push_back(end_id / 2);
	}
	if (imnodes::IsLinkDestroyed(&start_id)) {
		links_deleted.push_back(start_id);
	}
	if (imnodes::IsLinkHovered(&hover_id)) {
		links_hovered.push_back(hover_id);
		ImGui::OpenPopup("Link");
	}
	
	if (ImGui::BeginPopup("Link")) {
		if (links_hovered.empty()) ImGui::CloseCurrentPopup();
		new_weight.reset();
		if (
			!network.connections.empty() &&
			ImGui::SliderFloat("Weight", &network.connections[hover_id].w, 0, 1)
		) {
			new_weight = network.connections[hover_id].w;
		}

		ImGui::EndPopup();
	}
}


void Genome_Window::embed_render(Genome& genome) noexcept {
	ImGui::Columns(3);

	size_t n_inputs = genome.input_nodes;
	size_t n_outputs = genome.output_nodes;
	size_t n_hidden = genome.node_genes.size() - n_inputs - n_outputs;

	ImGui::Text("%zu inputs", n_inputs);
	ImGui::SameLine();
	if (ImGui::Button("+##in")) genome.add_new_input();

	ImGui::NextColumn();

	ImGui::Text("%zu outputs", n_outputs);
	ImGui::SameLine();
	if (ImGui::Button("+##out")) genome.add_new_output();

	ImGui::NextColumn();

	ImGui::Text("%zu hidden", n_hidden);
	ImGui::SameLine();
	if (ImGui::Button("+##hid")) genome.add_new_hidden();

	ImGui::Columns(1);

	ImGui::Text("Node genes");
	ImGui::BeginChildFrame(2, { ImGui::GetContentRegionAvail().x, 40 });

	for (auto& x : genome.node_genes) {
		ImGui::BeginGroup();
		ImGui::Text("Node %d", x.node_id);
		ImGui::Text(
			"%s", x.kind == x.Sensor ? "Input" : (x.kind == x.Output ? "Output" : "Hidden")
		);
		ImGui::EndGroup();
		ImGui::SameLine();
	}
	ImGui::Dummy({});
	ImGui::EndChildFrame();

	ImGui::Text("Connection genes");
	ImGui::BeginChildFrame(3, { ImGui::GetContentRegionAvail().x, 60 });

	for (auto& x : genome.connect_genes) {
		ImVec4 color;
		if (x.enabled) color = { 255, 255, 255, 255 };
		if (!x.enabled) color = { 100, 100, 100, 100 };

		ImGui::BeginGroup();
		ImGui::TextColored(color, "%d -> %d", x.in_id, x.out_id);
		ImGui::TextColored(color, "%4.2f", x.w);
		ImGui::TextColored(color, "Innov %d", x.innov);
		ImGui::EndGroup();
		ImGui::SameLine();
	}
	ImGui::Dummy({});

	ImGui::EndChildFrame();
}

void update_genome(Genome& gen, const Network_Window& network_window_input) noexcept {
	for (size_t i = 0; i < network_window_input.links_created.size(); i += 2) {
		gen.add_connection(
			network_window_input.links_created[i], network_window_input.links_created[i + 1], 0
		);
	}

	for (auto& x : network_window_input.links_deleted) {
		gen.del_connection(x);
	}

	if (network_window_input.new_weight) for (auto& x : network_window_input.links_hovered) {
		gen.update_weight(x, *network_window_input.new_weight);
	}
}

void Neat_Window::render(Neat& neat) noexcept {
	ImGui::Begin("Neat");


	if (ImGui::Button("Edit initial network")) {
		open_edit_initial_network = true;
	}

	ImGui::SameLine();

	if (ImGui::Button("Open Best Network")) {
		open_best_network = true;
	}
	if (open_best_network && ImGui::Begin("Best Network", &open_best_network)) {

		Genome* best_genome = nullptr;
		float best_fitness = 0;
		for (size_t i = 0; i < neat.population.size(); ++i) {
			if (neat.results[i].fitness > best_fitness) {
				best_fitness = neat.results[i].fitness;
				best_genome = &neat.population[i];
			}
		}

		if (best_genome) {
			auto net = best_genome->phenotype();
			best_network_window.embed_render(net);
		} else {
			ImGui::Text("There is no best genome.");
		}

		ImGui::End();
	}

	if (open_edit_initial_network && ImGui::Begin("Edit Network", &open_edit_initial_network)) {
		initial_genome_window.embed_render(initial_genome);
		ImGui::Separator();
		initial_network = initial_genome.phenotype();
		initial_network_window.embed_render(initial_network);
		update_genome(initial_genome, initial_network_window);
		ImGui::End();
	}
	ImGui::SameLine();
	if (ImGui::Button("Reset Initial Network")) {
		initial_genome = {};
		initial_network = {};
	}
	ImGui::Separator();

	if (ImGui::CollapsingHeader("Parameters")) {
		const char* eval_items[] = { "XOR" };
		int x = (int)evaluation;
		ImGui::ListBox("Eval", &x, eval_items, IM_ARRAYSIZE(eval_items));
		evaluation = (Evaluation)x;
		ImGui::SliderFloat("Kill Rate    ", &neat.survival_rate, 0, 1);
		ImGui::SliderFloat("Mutation Rate", &neat.mutation_rate, 0, 1);
		x = neat.population_size;
		ImGui::SliderInt("Population", &x, 0, 10'000);
		neat.population_size = x;
	}



	run_generation = ImGui::Button("Run a generation");
	if (ImGui::Button(auto_run ? "Stop" : "Run")) auto_run = !auto_run;


	size_t offset = 0;
	size_t size = max_fitness.size();
	if (max_fitness.size() > max_fitness_n_samples) {
		offset = max_fitness.size() - max_fitness_n_samples;
		size = max_fitness_n_samples;
	}

	ImGui::PlotLines(
		"Max Fitness",
		max_fitness.data() + offset,
		(int)size,
		0,
		nullptr,
		FLT_MAX,
		FLT_MAX,
		{500, 200}
	);

	int x = max_fitness_n_samples;
	ImGui::SliderInt("# Sample", &x, 0, 10000);
	max_fitness_n_samples = x;

	ImGui::End();
}


