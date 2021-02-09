#pragma once
#include "AI/Neat.hpp"
#include "Experiment.hpp"

extern Genome xor_genome() noexcept;
extern float xor_fitness(Network& net, void*) noexcept;
static Experiment xor_exp { "XOR", xor_genome, xor_fitness };

extern Genome spv_genome() noexcept;
extern void spv_render(Network& net, void* user) noexcept;
extern float spv_fitness(Network& net, void* user) noexcept;
static Experiment spv_exp { "SPV", spv_genome, spv_fitness, spv_render };

extern Genome dpv_genome() noexcept;
extern void dpv_render(Network& net, void* user) noexcept;
extern float dpv_fitness(Network& net, void* user) noexcept;
static Experiment dpv_exp { "DPV", dpv_genome, dpv_fitness, dpv_render };

extern Genome dp_genome() noexcept;
extern void dp_render(Network& net, void* user) noexcept;
extern float dp_fitness(Network& net, void* user) noexcept;
static Experiment dp_exp { "DP ", dp_genome, dp_fitness, dp_render };

extern Genome hdpv_genome() noexcept;
extern void hdpv_render(Network& net, void* user) noexcept;
extern float hdpv_fitness(Network& net, void* user) noexcept;
static Experiment hdp_exp { "HDPV", hdpv_genome, hdpv_fitness, hdpv_render };
