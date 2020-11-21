#include "Window.hpp"

#include "imgui.h"
#include "imgui/imnodes.h"

void Network_Window::embed_render(Network& network) noexcept {
	thread_local std::vector<bool> is_input_output;

	is_input_output.clear();
	is_input_output.resize(network.nodes.size());
	links_created.clear();
	imnodes::BeginNodeEditor();

	for (auto& x : network.input_nodes) {
		imnodes::BeginNode(x);

		imnodes::BeginNodeTitleBar();
		ImGui::Text("%s", "Input");
		imnodes::EndNodeTitleBar();

		ImGui::Text("%zu", x);

		imnodes::BeginOutputAttribute(network.nodes.size() + x);
		imnodes::EndOutputAttribute();

		imnodes::EndNode();

		is_input_output[x] = true;
	}

	for (auto& x : network.output_nodes) {
		imnodes::BeginNode(x);

		imnodes::BeginNodeTitleBar();
		ImGui::Text("%s", "Output");
		imnodes::EndNodeTitleBar();

		ImGui::Text("%zu", x);

		imnodes::BeginInputAttribute(x);
		imnodes::EndOutputAttribute();

		imnodes::EndNode();

		is_input_output[x] = true;
	}

	for (size_t i = 0; i < network.nodes.size(); ++i) if (!is_input_output[i]) {
		auto& x = network.nodes[i];

		imnodes::BeginNode(i);
		imnodes::BeginNodeTitleBar();
		ImGui::Text("%s", "Hidden");
		imnodes::EndNodeTitleBar();

		ImGui::Text("%zu", i);

		imnodes::BeginInputAttribute(i);
		imnodes::EndInputAttribute();

		imnodes::BeginOutputAttribute(network.nodes.size() + i);
		imnodes::EndOutputAttribute();

		imnodes::EndNode();
	}

	for (size_t i = 0; i < network.connections.size(); ++i) {
		imnodes::Link(
			i,
			network.connections[i].in + network.nodes.size(),
			network.connections[i].out
		);
	}

	imnodes::EndNodeEditor();

	int start_id = 0;
	int end_id = 0;
	if (imnodes::IsLinkCreated(&start_id, &end_id)) {
		links_created.push_back(start_id % network.nodes.size());
		links_created.push_back(end_id % network.nodes.size());
	}
	if (imnodes::IsLinkDestroyed(&start_id)) {
		links_deleted.push_back(start_id);
	}
}


void Genome_Window::embed_render(Genome& genome) noexcept {
	ImGui::Columns(2);

	size_t n_inputs = 0;
	size_t n_outputs = 0;

	for (auto& x : genome.node_genes) {
		if (x.kind == Genome::Node_Gene::Kind::Sensor) n_inputs++;
		if (x.kind == Genome::Node_Gene::Kind::Output) n_outputs++;
	}

	ImGui::Text("%zu inputs", n_inputs);
	ImGui::SameLine();
	if (ImGui::Button("+##in")) {
		genome.add_new_input();
	}
		

	ImGui::NextColumn();

	ImGui::Text("%zu outputs", n_outputs);
	ImGui::SameLine();
	if (ImGui::Button("+##out")) genome.add_new_output();

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
	for (size_t i = 0; i < network_window_input.links_created.size(); ++i) {
		gen.add_connection(
			network_window_input.links_created[i], network_window_input.links_created[i + 1], 0
		);
	}

	for (size_t i = 0; i < network_window_input.links_deleted.size(); ++i) {
		gen.del_connection(i);
	}
}

void render(Neat& neat) noexcept {
	ImGui::Begin("Neat");

	int x = neat.population.size();
	if (ImGui::SliderInt("Population", &x, 0, 1'000'000)) {
		if (x > neat.population.size()) x - neat.population.size();
		else                            neat.population.resize(x);
	}

	ImGui::Button("Run a generation");

	ImGui::End();
}

void Neat::fill_with(Genome g) noexcept {
	n_inputs = 0;
	n_outputs = 0;
	for (auto& x : g.node_genes) if (x.kind == x.Sensor) n_inputs++;
	for (auto& x : g.node_genes) if (x.kind == x.Output) n_outputs++;

	for (auto& x : population) x = g;
}

// >TODO(Tackwin):
void add_agents(size_t n) noexcept {}


