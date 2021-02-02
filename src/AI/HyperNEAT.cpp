#include "HyperNEAT.hpp"

void mark_node_plan(Substrate& s, size_t x, size_t y, Network::Node_Kind kind) noexcept {
	s.marked_idx[x + y * s.density] = kind;
}


Network hyper(std::vector<Substrate> substrates, Network& net) noexcept {
	Network hyper;
	std::unordered_map<
		size_t, std::unordered_map<size_t, float>
	> connections;
	std::unordered_map<size_t, Network::Node_Kind> kinds;

	size_t offset = 0;
	for (auto& x : substrates) {
		for (size_t x1 = 0; x1 < x.density; ++x1)
		for (size_t y1 = 0; y1 < x.density; ++y1) {
			auto id = offset + x1 + y1 * x.density;
			kinds[id] = Network::Node_Kind::Hidden;
			if (x.marked_idx.count(x1 + y1 * x.density)) kinds[id] = x.marked_idx[x1 + y1 * x.density];
		}

		offset += x.density * x.density;
	}

	offset = 0;
	for (size_t i = 0; i + 1 < substrates.size(); ++i) {
		auto& curr = substrates[i + 0];
		auto& next = substrates[i + 1];

		for (size_t x1 = 0; x1 < curr.density; ++x1)
		for (size_t y1 = 0; y1 < curr.density; ++y1)
		
		for (size_t x2 = 0; x2 < next.density; ++x2)
		for (size_t y2 = 0; y2 < next.density; ++y2) {

			auto id1 = x1 + y1 * curr.density + offset;
			auto id2 = x2 + y2 * next.density + offset + curr.density * curr.density;

			float ins[] = {
				(float)x1 / curr.density,
				(float)y1 / curr.density,
				(float)x2 / next.density,
				(float)y2 / next.density,
				1.f
			};

			auto& v = connections[id1][id2];
			net.compute_clear(ins, 5, &v, 1);
		}

		offset += curr.density * curr.density;
	}
	if (substrates.size() > 0) offset += substrates.back().density * substrates.back().density;

	std::unordered_map<size_t, size_t> map_idx;
	for (auto& [x, k] : kinds) {
		if (k == Network::Node_Kind::Sensor) {
			map_idx[x] = hyper.input_nodes++;
		}
	}
	for (auto& [x, k] : kinds) {
		if (k == Network::Node_Kind::Output) {
			map_idx[x] = hyper.input_nodes + hyper.output_nodes++;
		}
	}
	size_t n_hidden = 0;
	for (auto& [x, k] : kinds) {
		if (k == Network::Node_Kind::Hidden) {
			map_idx[x] = hyper.input_nodes + hyper.output_nodes + n_hidden++;
		}
	}

	hyper.nodes.resize(offset);
	hyper.node_kinds.resize(offset, Network::Node_Kind::Hidden);

	for (auto& [a, b] : map_idx) hyper.node_kinds[b] = kinds[a];

	for (auto& [a, c] : connections) for (auto& [b, w] : c) {
		Network::Connections conn;
		conn.in  = (std::uint32_t)map_idx[a];
		conn.out = (std::uint32_t)map_idx[b];
		conn.w   = w;

		hyper.connections.push_back(conn);
	}

	return hyper;
}
