#pragma once

#include "Neat.hpp"

struct Substrate {
	enum Kind {
		Grid = 0,
		Count
	};

	Kind kind = Kind::Grid;
	size_t density = 5;

	std::unordered_map<size_t, Network::Node_Kind> marked_idx;
};

extern void mark_node_plan(Substrate& s, size_t x, size_t y, Network::Node_Kind kind) noexcept;

extern Network hyper(std::vector<Substrate> substrates, Network& Network) noexcept;
