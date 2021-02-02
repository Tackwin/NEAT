#pragma once
#include "AI/Neat.hpp"

extern float xor_fitness(Network& net, void*) noexcept;

extern void spv_render(Network& net, void* user) noexcept;
extern float spv_fitness(Network& net, void* user) noexcept;

extern void dpv_render(Network& net, void* user) noexcept;
extern float dpv_fitness(Network& net, void* user) noexcept;

extern void dp_render(Network& net, void* user) noexcept;
extern float dp_fitness(Network& net, void* user) noexcept;

extern void hdpv_render(Network& net, void* user) noexcept;
extern float hdpv_fitness(Network& net, void* user) noexcept;
