#include "Window.hpp"

#include "imgui.h"
#include "imgui/imnodes.h"

#include <algorithm>

void Network_Window::embed_render(Network& network) noexcept {
	links_created.clear();
	links_deleted.clear();
	links_hovered.clear();

	size_t n_inputs = 0;
	size_t n_outputs = 0;

	if (ImGui::Button("Reset")) network = {};
	ImGui::SameLine();
	if (ImGui::Button("Fully connect")) {
		for (size_t i = 0; i < network.nodes.size(); ++i)
			if (network.node_kinds[i] != Network::Node_Kind::Output)
			for (size_t j = i; j < network.nodes.size(); ++j)
				if (network.node_kinds[j] != Network::Node_Kind::Sensor) {
					Network::Connections c;
					c.in = i;
					c.out = j;
					c.w = 1;
					bool not_found = true;
					for (auto& x : network.connections) if (x.in == c.in && x.out == c.out) {
						not_found = false;
					}

					if (not_found) {
						links_created.push_back(c.in);
						links_created.push_back(c.out);
					}
				}

	}

	imnodes::BeginNodeEditor();

	for (size_t i = 0; i < network.nodes.size(); ++i) {
		if (network.node_kinds[i] == Network::Node_Kind::Sensor) n_inputs++;
		if (network.node_kinds[i] == Network::Node_Kind::Output) n_outputs++;
		
		imnodes::BeginNode(i);

		imnodes::BeginNodeTitleBar();
		ImGui::Text(
			"%s %zu",
			network.node_kinds[i] == Network::Node_Kind::Hidden ? "Hidden" :
			(network.node_kinds[i] == Network::Node_Kind::Output ? "Output" :
			"Input"),
			i
		);
		imnodes::EndNodeTitleBar();

		ImGui::Text("%s", node_activation_func_str(network.nodes[i].f));

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


	if (ImGui::CollapsingHeader("Test")) {
		thread_local std::vector<float> inputs;
		thread_local std::vector<float> outputs;

		inputs.resize(n_inputs);
		outputs.resize(n_outputs);

		ImGui::Columns(2);

		for (auto& x : inputs) {
			ImGui::PushID(&x);
			ImGui::SliderFloat("Input", &x, 0, 1);
			ImGui::PopID();
		}

		ImGui::NextColumn();

		network.compute(inputs.data(), inputs.size(), outputs.data(), outputs.size());

		for (auto& x : outputs) {
			ImGui::PushID(&x);
			ImGui::Text("Output %f", x);
			ImGui::PopID();
		}

		ImGui::Columns(1);
	}
}


void Genome_Window::embed_render(Genome& genome) noexcept {
	if (ImGui::Button("Save")) {
		auto data = genome.serialize();

		FILE* f = fopen("init_genome", "wb");
		fwrite(data.data(), 1, data.size(), f);
		fclose(f);
	}
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
			network_window_input.links_created[i], network_window_input.links_created[i + 1], 1
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

	if (ImGui::Button("Reset NEAT")) {
		neat = {};
	}

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
		x = neat.population_size;
		ImGui::SliderInt("Population", &x, 0, 10'000);
		neat.population_size = x;

		if (ImGui::CollapsingHeader("Rates")) {
			ImGui::SliderFloat("Kill Rate    ", &neat.survival_rate, 0, 1);
			ImGui::SliderFloat("Mutation Rate", &neat.mutation_rate, 0, 1);
			ImGui::SliderFloat("Young advantage", &neat.young_advantage, 0, 1);
			ImGui::SliderFloat("Specie Crowd Rate", &neat.specie_crowd_rate, 0, 1);
			ImGui::SliderFloat("Complexity cost", &neat.complexity_cost, 0, 1);
		}

		if (ImGui::CollapsingHeader("Specie Diff Parameters")) {
			ImGui::SliderFloat("C1", &neat.c_1, 0, 1);
			ImGui::SliderFloat("C2", &neat.c_2, 0, 1);
			ImGui::SliderFloat("C3", &neat.c_3, 0, 1);
			if (neat.specifie_number_of_species) {
				x = neat.preferred_number_of_species;
				ImGui::SliderInt("#Specie", &x, 0, neat.population_size);
				neat.preferred_number_of_species = x;
			} else {
				ImGui::SliderFloat("Dt", &neat.specie_dt, 0, 10);
			}
			ImGui::SameLine();
			ImGui::Checkbox("Hard", &neat.specifie_number_of_species);
		}
	}

	run_generation = ImGui::Button("Run a generation");
	if (ImGui::Button(auto_run ? "Stop" : "Run")) auto_run = !auto_run;

	if (ImGui::CollapsingHeader("Infos")) {
		ImGui::Text("Generation: % 10lld\n", (long long int)neat.generation_number);
	}

	size_t offset = 0;
	size_t size = max_fitness.size();
	if (max_fitness.size() > max_fitness_n_samples) {
		offset = max_fitness.size() - max_fitness_n_samples;
		size = max_fitness_n_samples;
	}

	thread_local bool adjusted = false;
	ImGui::Checkbox("Adjusted", &adjusted);

	auto to_plot = max_fitness.data();
	if (adjusted) to_plot = max_adjusted_fitness.data();
	ImGui::PlotLines(
		"Max Fitness",
		to_plot + offset,
		(int)size,
		0,
		nullptr,
		0,
		FLT_MAX,
		{500, 200}
	);
	if (fitness_histograms.size() > 0) {
		to_plot = fitness_histograms.back().data();
		if (adjusted) to_plot = adjusted_fitness_histograms.back().data();

		ImGui::PlotHistogram("Fitness dist", to_plot, 100, 0, nullptr, 0, FLT_MAX, {500, 200});
	}
	ImGui::PlotHistogram(
		"Species",
		[](void* data, int idx) {
			return (float)((Neat::Specie*)data)[idx].size;
		},
		neat.species.data(),
		neat.species.size(),
		0,
		nullptr,
		0,
		FLT_MAX,
		{500, 200}
	);

	int x = max_fitness_n_samples;
	ImGui::SliderInt("# Sample", &x, 0, 10000);
	max_fitness_n_samples = x;

	ImGui::End();
}


void Neat_Window::get_stats(const std::vector<Neat::Result>& results) noexcept {
	thread_local std::array<float, 100> hist;
	auto best_adjusted_fitness = 0.f;
	auto best_fitness = 0.f;

	for (auto& x : results) {
		best_adjusted_fitness = std::max(best_adjusted_fitness, x.adjusted_fitness);
		best_fitness = std::max(best_fitness, x.fitness);
	}

	memset(hist.data(), 0, hist.size());
	for (auto& x : results) {
		auto t = (size_t)(99 * x.fitness / best_fitness);

		hist[t]++;
	}
	max_fitness.push_back((float)best_fitness);
	fitness_histograms.push_back(hist);

	memset(hist.data(), 0, hist.size());
	for (auto& x : results) {
		auto t = (size_t)(99 * x.adjusted_fitness / best_adjusted_fitness);

		hist[t]++;
	}
	max_adjusted_fitness.push_back((float)best_adjusted_fitness);
	adjusted_fitness_histograms.push_back(hist);
}
