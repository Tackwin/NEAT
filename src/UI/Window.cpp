#include "Window.hpp"

#define IMGUI_DEFINE_MATH_OPERATORS
#include "imgui.h"
#include "imgui_internal.h"
#include "imgui/imnodes.h"

#include "xstd.hpp"

#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <stdio.h>

namespace ImGui {
	bool SliderSize(
		const char* label,
		size_t* v,
		size_t v_min,
		size_t v_max,
		const char* format = "%d",
		ImGuiSliderFlags flags = 0
	) {
		int x = *v;
		bool ret = SliderInt(label, &x, (int)v_min, (int)v_max, format, flags);
		*v = x;
		return ret;
	}
}


void Network_Window::embed_render(Network& network) noexcept {
	links_created.clear();
	links_deleted.clear();
	links_hovered.clear();

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
	ImGui::SameLine();
	layout = ImGui::Button("Layout");

	imnodes::EditorContextSet(ctx.ctx);
	
	imnodes::BeginNodeEditor();

	for (size_t i = 0; i < network.nodes.size(); ++i) {
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

		imnodes::BeginOutputAttribute(i * 2 + 0);
		imnodes::EndOutputAttribute();
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
	hovered = imnodes::IsEditorHovered();

	if (layout) auto_layout(network);

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

		inputs.resize(network.input_nodes);
		outputs.resize(network.output_nodes);

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


void Network_Window::auto_layout(Network& network) noexcept {
	// Honestly i could search for the right layout with a good cost function. Maybe i should...

	struct Vec2f {
		float x = 0;
		float y = 0;
	};

	std::vector<Vec2f> positions;
	positions.resize(network.nodes.size());

	std::unordered_map<size_t, std::vector<size_t>> connections;
	for (auto& x : network.connections) connections[x.in].push_back(x.out);

	std::vector<std::pair<size_t, size_t>> open;
	std::unordered_set<size_t> closed;

	std::vector<size_t> layers;
	layers.resize(network.nodes.size());

	size_t n_layer = 0;

	for (size_t i = 0; i < network.input_nodes; ++i) open.push_back({i, 0});

	while (!open.empty()) {
		auto [it, l] = open.back();
		open.pop_back();

		if (closed.find(it) != std::end(closed)) continue;
		closed.insert(it);

		n_layer = std::max(n_layer, l);
		layers[it] = l;

		for (auto& x : connections[it]) open.push_back({x, l+1});
	}

	std::vector<size_t> n_by_layers;
	n_by_layers.resize(n_layer + 1);
	for (auto& x : layers) n_by_layers[x]++;

	std::vector<size_t> max_by_layers = n_by_layers;

	auto max_width = std::max_element(std::begin(n_by_layers), std::end(n_by_layers));

	for (size_t i = 0; i < network.nodes.size(); ++i) {
		auto l = layers[i];

		imnodes::SetNodeGridSpacePos(i, {
			100.f * l,
			100.f * n_by_layers[l] - 50.f * max_by_layers[l]
		});

		n_by_layers[layers[i]]--;
	}

}


void Genome_Window::embed_render(Genome& genome) noexcept {
	if (ImGui::Button("Save")) {
		auto data = genome.serialize();

		FILE* f = fopen("init_genome", "wb");
		fwrite(data.data(), 1, data.size(), f);
		fclose(f);
	}
	ImGui::SameLine();
	if (ImGui::Button("Load")) {
		FILE* f = fopen("init_genome", "rb");

		fseek(f, 0, SEEK_END);
		long fsize = ftell(f);
		fseek(f, 0, SEEK_SET);  /* same as rewind(f); */

		std::vector<std::uint8_t> data;
		data.resize(fsize);
		fread(data.data(), 1, fsize, f);
		fclose(f);

		genome = Genome::deserialize(data);
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
	ImGui::SameLine();
	if (ImGui::Button("Open Best Phenotype")) {
		open_best_phenotype = true;
	}

	if (open_best_network) {
		if (ImGui::Begin("Best Network", &open_best_network)) {
			auto net = best_genome.phenotype();

			best_genome_window.embed_render(best_genome);
			ImGui::Separator();
			best_network_window.embed_render(net);
		}
		ImGui::End();
	}

	if (open_best_phenotype) render_best(neat);

	if (open_edit_initial_network) {
		if (ImGui::Begin("Edit Network", &open_edit_initial_network)) {
			initial_genome_window.embed_render(initial_genome);
			ImGui::Separator();
			initial_network = initial_genome.phenotype();
			initial_network_window.embed_render(initial_network);
			update_genome(initial_genome, initial_network_window);
		}
		ImGui::End();
	}
	ImGui::SameLine();
	if (ImGui::Button("Reset Initial Network")) {
		initial_genome = {};
		initial_network = {};
	}
	ImGui::Separator();

	if (ImGui::CollapsingHeader("Parameters")) {
		const char* eval_items[] = {
			"XOR",
			"Single Pole Velocity",
			"Double Pole Velocity",
			"Double Pole",
			"Hyper"
		};
		int x = (int)evaluation;
		ImGui::ListBox("Eval", &x, eval_items, 5);
		evaluation = (Evaluation)x;
		ImGui::SliderSize("Population", &neat.population_size, 0, 100'000);

		float indent = 64.f;
		ImGui::PushItemWidth(200);
		ImGui::Indent(indent);
		if (ImGui::CollapsingHeader("Rates")) {
			ImGui::Indent(indent);
			#define X(x) ImGui::SliderFloat(#x, &x, 0, 1);
			X(neat.survival_rate);
			X(neat.young_advantage);
			X(neat.aging_specie_penality);
			X(neat.specie_crowd_rate);
			X(neat.complexity_cost);
			X(neat.population_competition_rate);
			X(neat.mutation_rate);
			ImGui::SliderSize("Sterilize specie", &neat.age_cutoff_sterile_specie, 0, 1000);
			if (ImGui::CollapsingHeader("Mutation")) {
				ImGui::Indent(indent);
				X(neat.params.update_connection_rate);
				X(neat.params.toggle_connection_rate);
				X(neat.params.new_connection_weight_rate);
				X(neat.params.add_connection_rate);
				X(neat.params.update_node_act_rate);
				X(neat.params.add_node_rate);
				X(neat.params.del_node_rate);
				ImGui::Unindent(indent);
			#undef X
			}
			ImGui::Unindent(indent);
		}
		ImGui::PopItemWidth();

		if (ImGui::CollapsingHeader("Specie Diff Parameters")) {
			ImGui::Indent(indent);
			ImGui::SliderFloat("C1", &neat.c_1, 0, 1);
			ImGui::SliderFloat("C2", &neat.c_2, 0, 1);
			ImGui::SliderFloat("C3", &neat.c_3, 0, 1);
			ImGui::SliderSize("Specie size", &neat.min_specie_size_advantage, 0, 100);
			if (neat.specifie_number_of_species) {
				ImGui::SliderSize(
					"#Specie", &neat.preferred_number_of_species, 0, neat.population_size
				);
			} else {
				ImGui::SliderFloat("Dt", &neat.specie_dt, 0, 10);
			}
			ImGui::SameLine();
			ImGui::Checkbox("Hard", &neat.specifie_number_of_species);
			ImGui::Unindent(indent);
		}
		ImGui::Unindent(indent);
	}
	if (ImGui::CollapsingHeader("Visualisation")) {
		ImGui::Checkbox("Capture population", &capture_population);
		if (ImGui::Button("Explore population")) open_explore_population = true;
	}
	if (ImGui::CollapsingHeader("Performance")) {
		ImGui::SliderSize("Threads", &neat.threads_to_use, 1, 12);
	}

	if (open_explore_population) population_window.render(population_snapshots, results_snapshots);

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

	float baseline = 0;
	for (size_t i = 0; i < size; ++i) baseline = std::min(to_plot[offset + i], baseline);

	ImGui::PlotLines(
		"Max Fitness",
		to_plot + offset,
		(int)size,
		0,
		nullptr,
		baseline,
		FLT_MAX,
		{500, 200}
	);
	ImGui::SameLine();
	ImGui::Text("%f", max_fitness.empty() ? 0 : max_fitness.back());
	if (fitness_histograms.size() > 0) {
		to_plot = fitness_histograms.back().data();
		if (adjusted) to_plot = adjusted_fitness_histograms.back().data();

		ImGui::PlotHistogram("Fitness dist", to_plot, 100, 0, nullptr, 0, FLT_MAX, {500, 200});
	}

	thread_local std::vector<size_t> xs;
	xs.clear();
	for (size_t i = 0; i < neat.generation_number; ++i) xs.push_back(i);

	if (!species_infos.empty()) {
		open_species_visu |= ImGui::Button("Open species visu");
		if (open_species_visu) visu_species_size();
	}

	ImGui::SliderSize("# Sample", &max_fitness_n_samples, 0, 10000);

	ImGui::End();
}

void Neat_Window::render_best(Neat& neat) noexcept {
	ImGui::Begin("Best Phenotype", &open_best_phenotype);
	if (render_phenotype) {
		auto n = best_genome.phenotype();
		render_phenotype(n, nullptr);
	} else {
		ImGui::Text("No visualisation is available for this problem. :(");
	}
	ImGui::End();
}

void Neat_Window::get_stats(const Neat& neat, const std::vector<Neat::Result>& results) noexcept {
	thread_local std::array<float, 100> hist;
	auto worst_adjusted_fitness = FLT_MAX;
	auto worst_fitness = FLT_MAX;
	auto best_adjusted_fitness = -FLT_MAX;
	auto best_fitness = -FLT_MAX;

	for (auto& x : results) {
		if (best_fitness < x.fitness) {
			best_fitness = x.fitness;
			best_genome = *x.g;
		}

		best_adjusted_fitness = std::max(best_adjusted_fitness, x.adjusted_fitness);
		best_fitness = std::max(best_fitness, x.fitness);
		worst_adjusted_fitness = std::min(worst_adjusted_fitness, x.adjusted_fitness);
		worst_fitness = std::min(worst_fitness, x.fitness);
	}

	for (auto& x : hist) x = 0;
	for (auto& x : results) {
		size_t t = 99;
		if (x.fitness != best_fitness)
			t = (size_t)(99 * (x.fitness - worst_fitness) / (best_fitness - worst_fitness));

		hist[t]++;
	}
	max_fitness.push_back((float)best_fitness);
	fitness_histograms.push_back(hist);

	for (auto& x : hist) x = 0;
	for (auto& x : results) {
		size_t t = 99;
		if (x.adjusted_fitness != best_adjusted_fitness)
			t = (size_t)(99 *
				(x.adjusted_fitness - worst_adjusted_fitness) /
				(best_adjusted_fitness - worst_adjusted_fitness)
			);

		hist[t]++;
	}
	max_adjusted_fitness.push_back((float)best_adjusted_fitness);
	adjusted_fitness_histograms.push_back(hist);

	if (capture_population) {
		population_snapshots.push_back(neat.population);
		results_snapshots.push_back(neat.results);
	}

	std::vector<Specie_Info> infos;
	for (auto& x : neat.species) {
		Specie_Info i;
		i.n = x.size;
		i.id = x.num;
		infos.push_back(i);
	}
	species_infos.push_back(infos);
}


void Population_Window::render(
	const std::vector<std::vector<Genome>>& pop,
	const std::vector<std::vector<Neat::Result>>& res
) noexcept {
	if (ImGui::Begin("Population"))
		embed_render(pop, res);
	ImGui::End();
}

void Population_Window::embed_render(
	const std::vector<std::vector<Genome>>& pop,
	const std::vector<std::vector<Neat::Result>>& res
) noexcept {
	if (pop.empty()) {
		ImGui::Text("Start running generations.");
		return;
	}
	ImGui::SliderSize("Gen", &generation, 0, pop.size() - 1);
	generation = std::max((size_t)0, std::min(generation, pop.size() - 1));
	
	float min_fitness_possible = FLT_MAX;
	float max_fitness_possible = -FLT_MAX;
	size_t min_neurons_possible = SIZE_MAX;
	size_t max_neurons_possible = 0;
	for (auto& x : res[generation]) {
		min_fitness_possible = std::min(min_fitness_possible, x.fitness);
		max_fitness_possible = std::max(max_fitness_possible, x.fitness);
	}
	for (auto& x : pop[generation]) {
		min_neurons_possible = std::min(min_neurons_possible, x.node_genes.size());
		max_neurons_possible = std::max(max_neurons_possible, x.node_genes.size());
	}

	ImGui::PushItemWidth(100);
	ImGui::SliderFloat("Min Ftiness", &min_fitness, min_fitness_possible, max_fitness);
	ImGui::SameLine();
	ImGui::SliderFloat("Max Fitness", &max_fitness, min_fitness, max_fitness_possible);

	ImGui::SliderSize("Min Neurons", &min_neurons, min_neurons_possible, max_neurons);
	ImGui::SameLine();
	ImGui::SliderSize("Max Neurons", &max_neurons, min_neurons, max_neurons_possible);
	ImGui::PopItemWidth();

	min_neurons = std::min(max_neurons_possible, std::max(min_neurons_possible, min_neurons));
	max_neurons = std::min(max_neurons_possible, std::max(min_neurons_possible, max_neurons));
	min_fitness = std::min(max_fitness_possible, std::max(min_fitness_possible, min_fitness));
	max_fitness = std::min(max_fitness_possible, std::max(min_fitness_possible, max_fitness));

	ImGui::BeginChildFrame(1, ImGui::GetContentRegionAvail());

	network_windows.resize(pop[generation].size());


	size_t n = ImGui::GetContentRegionAvail().x / 150;

	thread_local std::unordered_map<size_t, ImVec2> net_size;
	thread_local std::unordered_map<size_t, bool> net_focused;

	for (size_t i = 0; i < pop[generation].size(); ++i) {
		auto& it = pop[generation][i];
		auto& it_res = res[generation][i];
		if (!(min_fitness <= it_res.fitness && it_res.fitness <= max_fitness)) continue;
		if (!(min_neurons <= it.node_genes.size() && it.node_genes.size() <= max_neurons)) continue;

		if (net_size.count(i) == 0) net_size[i] = { 125, 125 };
		if (net_focused[i]) net_size[i] = { 600, 600 };

		ImGui::BeginChild(i + 2, net_size[i], true);
		bool b = net_focused[i];
		ImGui::Checkbox("Focused", &b);
		net_focused[i] = b;

		ImGui::SameLine();

		ImGui::Text(
			"Fitness % 10.5f Adjusted fitness % 10.5f", it_res.fitness, it_res.adjusted_fitness
		);

		auto phen = pop[generation][i].phenotype();
		network_windows[i].embed_render(phen);

		if (ImGui::IsItemHovered(ImGuiHoveredFlags_RectOnly) || network_windows[i].hovered) net_size[i] = { 600, 600 };
		else                                                        net_size[i] = { 150, 150 };
		ImGui::EndChild();
		if (n > 0 && (i % n) != (n - 1)) ImGui::SameLine();
	}
	ImGui::Dummy({0, 0});

	ImGui::EndChildFrame();
}

void Neat_Window::visu_species_size() noexcept {
	if (!ImGui::Begin("Species", &open_species_visu)) {
		ImGui::End();
		return;
	}

	auto draw_list = ImGui::GetWindowDrawList();
	auto offset = ImGui::GetWindowPos();
	auto pos = ImGui::GetCursorPos();
	auto size = ImGui::GetContentRegionAvail();
	auto center = offset + size / 2;

	size_t cutoff_size = 1;

	if (species_infos.empty()) {
		ImGui::Text("Try running some generations.");
		ImGui::End();
		return;
	}

	auto& infos = species_infos.back();

	struct Circle {
		size_t specie_idx = 0;
		ImVec2 pos        = {};
		ImColor color     = {};
		float size        = .0f;
		float target_size = .0f;
	};
	thread_local std::unordered_map<size_t, Circle> last;
	std::unordered_map<size_t, Circle> circles;
	auto hash = [] (size_t x) { return (uint32_t)(x * 2654435761); };

	for (size_t i = 0; i < infos.size(); ++i) if (infos[i].n > cutoff_size) {
		auto& x = infos[i];
		Circle c;
		c.pos = { (float)(randomd() * 200 - 100), (float)(randomd() * 200 - 100)};
		if (last.find(x.id) != std::end(last)) c.pos = last[x.id].pos;
		if (last.find(x.id) != std::end(last)) c.size = last[x.id].size;
		c.color = hash(x.id);
		c.color.Value.w = 255;
		c.target_size = sqrtf((float)x.n);
		c.specie_idx = i;
		circles[x.id] = c;
	}

	float center_pression = 200;
	float outward_pression = 600;
	float growth_size = 4.f;
	auto dt = ImGui::GetIO().DeltaTime;
	dt = dt > 0.01f ? 0.01f : dt;

	for (auto& [_, x] : circles) {
		x.size += growth_size * dt * (x.target_size - x.size);
	}

	// we cetner every one.
	for (auto& [_, x] : circles) {
		float d = x.pos.x * x.pos.x + x.pos.y * x.pos.y;
		d = sqrtf((float)d) + 0.01f;

		auto pull = (d > center_pression) ? center_pression : d;

		x.pos.x -= dt * pull * x.pos.x / d;
		x.pos.y -= dt * pull * x.pos.y / d;
	}

	for (auto& [a_id, a] : circles) {
		for (auto& [b_id, b] : circles) if (a_id < b_id) {
			float d = (float)(
				(a.pos.x - b.pos.x) * (a.pos.x - b.pos.x) +
				(a.pos.y - b.pos.y) * (a.pos.y - b.pos.y)
			);
			d = sqrtf((float)d) + 0.01f;

			auto denom = d / (a.size + b.size);
			denom *= denom;
			denom *= denom;
			auto f = outward_pression / (1 + denom);

			a.pos.x += dt * f * (a.pos.x - b.pos.x) / (d + a.size);
			a.pos.y += dt * f * (a.pos.y - b.pos.y) / (d + a.size);
			b.pos.x += dt * f * (b.pos.x - a.pos.x) / (d + b.size);
			b.pos.y += dt * f * (b.pos.y - a.pos.y) / (d + b.size);
		}
	}
	last = circles;

	for (auto& [_, x] : circles) {
		draw_list->AddCircleFilled(center + x.pos, x.size, x.color);
	}

	auto mouse_pos = ImGui::GetMousePos();
	for (auto& [id, x] : circles) {
		auto d = mouse_pos - (center + x.pos);
		auto dist = sqrtf((float)(d.x * d.x + d.y * d.y));

		if (dist > x.size) continue;

		ImGui::BeginTooltip();
		ImGui::Text("Specie #: % 5d", (int)infos[x.specie_idx].id);
		ImGui::Text("Size    : % 5d", (int)infos[x.specie_idx].n);
		ImGui::EndTooltip();
	}

	ImGui::End();
}


