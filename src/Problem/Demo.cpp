#include "Demo.hpp"

#include "AI/HyperNEAT.hpp"

#define IMGUI_DEFINE_MATH_OPERATORS
#include "imgui/imgui.h"
#include "imgui/imgui_internal.h"

#include "xstd.hpp"

#include <algorithm>


float xor_fitness(Network& net, void*) noexcept {
	float inputs[] = {
		0, 0, 1,
		1, 0, 1,
		0, 1, 1,
		1, 1, 1
	};
	float target_outputs[] = {
		0, 1, 1, 0
	};
	float predicted_outputs[4];

	net.compute_clear(inputs + 0, 3, predicted_outputs + 0, 1);
	net.compute_clear(inputs + 3, 3, predicted_outputs + 1, 1);
	net.compute_clear(inputs + 6, 3, predicted_outputs + 2, 1);
	net.compute_clear(inputs + 9, 3, predicted_outputs + 3, 1);

	for (auto& x : predicted_outputs) x = std::max(0.f, std::min(1.f, x));

	float s = 0;
	for (size_t i = 0; i < 4; ++i) {
		auto dt = (predicted_outputs[i] - target_outputs[i]);
		s += dt * dt;
	}

	return 1 / (0.1f + s);
}


// ===========================  SPV
// ===========================  SPV
// ===========================  SPV
// ===========================  SPV
// ===========================  SPV


struct Single_Pole_State {
	double cart_x = 0;
	double cart_v = 0;

	double pole_θ = 0;
	double pole_v = 0;
};
auto sp_compute = [](Single_Pole_State state, Network& net) -> float {
	float inputs[] = {
		(float)state.cart_x,
		(float)state.cart_v,
		(float)state.pole_θ,
		(float)state.pole_v,
		1.f // biais
	};

	inputs[0] = (inputs[0] + 2.4f) / 4.8f;
	inputs[1] = sig(inputs[1]);
	inputs[2] = (atan2(sinf(inputs[2]), cosf(inputs[2])) + M_PI) / (2 * M_PI);
	inputs[3] = sig(inputs[3]);

	float outputs[1];

	net.compute(inputs, 5, outputs, 1);

	return outputs[0];
};

auto sp_timestep = [](Single_Pole_State in, double f) -> Single_Pole_State {
	Single_Pole_State out;
	constexpr double m    = 0.1;    // mass of the pole
	constexpr double m_c  = 1.0;    // mass of the cart
	constexpr double g    = 9.8;    // gravity
	constexpr double l    = 0.2;    // length of the pole
	constexpr double mu_c = 0.1;    // cart firction coeff
	constexpr double mu_p = .003;   // joint firction coeff
	constexpr double M    = m + m_c; // combined mass of the system
	constexpr double t_l  = 4.8;    // trach width in meter

	constexpr double dt = 0.01; // timestep in second

	// Extract state into named variables.
	// Precompute some reused values.
	double sin_theta = sin(in.pole_θ);
	double cos_theta = cos(in.pole_θ);
	double cos_theta_sqr = cos_theta * cos_theta;
	double thetav_sqr = in.pole_v * in.pole_v;

	double l_hat = l / 2;

	if (out.cart_x == +t_l / 2 && f > 0) f = 0;
	if (out.cart_x == -t_l / 2 && f < 0) f = 0;

	// Calc cart horizontal acceleration.
	auto xa = (
		m * g * sin_theta * cos_theta -
		(7.0 / 3.0) * (f + m * l_hat * thetav_sqr * sin_theta - mu_c * in.cart_v) -
		((mu_p * in.pole_v * cos_theta) / l_hat)
	) / (m * cos_theta_sqr - (7.0 / 3.0) * M);

	// Calc pole angular acceleration.
	auto θa =
		3.0 / (7.0 * l_hat) *
		(g * sin_theta - xa * cos_theta - (mu_p * in.pole_v) / (m * l_hat));

	out.cart_x = in.cart_x + in.cart_v * dt;
	out.cart_v = in.cart_v + xa * dt;
	out.pole_θ = in.pole_θ + in.pole_v * dt;
	out.pole_v = in.pole_v + θa * dt;

	if (out.cart_x < -t_l / 2) {
		out.cart_x = -t_l / 2;
		out.cart_v = 0;
	}
	
	if (out.cart_x > +t_l / 2) {
		out.cart_x = +t_l / 2;
		out.cart_v = 0;
	}

	return out;
};

void spv_render(Network& net, void* user) noexcept {
	//Single_Pole_State& state = *(Single_Pole_State*)user;
	std::vector<Single_Pole_State> results;
	Single_Pole_State state;
	state.cart_v = 0;
	state.cart_x = 0;
	state.pole_v = 0;
	state.pole_θ = (0.1f / 180) * M_PI;

	size_t N_T = 1000;

	results.push_back(state);

	for (size_t i = 0; i < N_T; i += 2) {
		auto f = sp_compute(results.back(), net) * 20 - 10;
		results.push_back(sp_timestep(results.back(), f));
		results.push_back(sp_timestep(results.back(), 0));
	}

	thread_local float time = 0;
	thread_local bool run = false;

	if (ImGui::Button(run ? "Pause" : " Run ")) run = !run;
	ImGui::SameLine();
	ImGui::SliderFloat("Time", &time, 0, 1);

	if (run) {
		time += ImGui::GetIO().DeltaTime / 10;
	}

	time = fmodf(time, 1.f);
	time = std::max(time, 0.f);

	state = results[time * (results.size() - 1)];

	auto draw_list = ImGui::GetWindowDrawList();

	auto offset = ImGui::GetWindowPos();
	auto pos = ImGui::GetCursorPos();
	auto size = ImGui::GetContentRegionAvail();
	auto center = pos + size / 2;

	auto left_end = ImVec2{ pos.x + size.x / 20, pos.y + 4 * size.y / 5 };
	auto right_end = ImVec2{ pos.x + 19 * size.x / 20, pos.y + 4 * size.y / 5 };


	auto scale = (right_end.x - left_end.x) / 4.8f;

	auto t = ((float)state.cart_x + 2.4f) / 4.8f;
	auto cart_pos = ImVec2{ left_end.x + (right_end.x - left_end.x) * t, left_end.y };

	auto pole_start = cart_pos;
	auto pole_end = pole_start + ImVec2{
		cosf(state.pole_θ - M_PI_2) * 0.2f * scale,
		sinf(state.pole_θ - M_PI_2) * 0.2f * scale
	};

	// draw left end
	draw_list->AddRectFilled(
		offset + left_end - ImVec2(0.05f * scale, 0.05f * scale),
		offset + left_end + ImVec2(0.05f * scale, 0.05f * scale),
		IM_COL32(255, 0, 0, 255)
	);
	// draw right end
	draw_list->AddRectFilled(
		offset + right_end - ImVec2(0.05f * scale, 0.05f * scale),
		offset + right_end + ImVec2(0.05f * scale, 0.05f * scale),
		IM_COL32(255, 0, 0, 255)
	);
	// draw base
	draw_list->AddRectFilled(
		offset + left_end  - ImVec2(0.0125f * scale, 0.0125f * scale),
		offset + right_end + ImVec2(0.0125f * scale, 0.0125f * scale),
		IM_COL32(255, 0, 0, 120)
	);
	// draw probe
	draw_list->AddRectFilled(
		offset + cart_pos + ImVec2(0.07f * scale, 0.07f * scale),
		offset + cart_pos - ImVec2(0.07f * scale, 0.07f * scale),
		IM_COL32(0, 255, 0, 255)
	);
	// draw pole
	draw_list->AddLine(
		offset + pole_start,
		offset + pole_end,
		IM_COL32(0, 0, 255, 200),
		0.02f * scale
	);

	ImGui::SetCursorPos(pos + size);
}

float spv_fitness(Network& net, void*) noexcept {
	thread_local std::vector<Single_Pole_State> results;
	results.clear();

	Single_Pole_State state;
	state.cart_v = 0;
	state.cart_x = 0;
	state.pole_v = 0;
	state.pole_θ = (0.1f / 180) * M_PI;

	size_t N_T = 1000;

	results.push_back(state);
	auto is_upright = [](Single_Pole_State x) {
		return std::abs(x.pole_θ) < M_PI_4 / 6;
	};

	size_t t = 0;
	for (; t < N_T; t += 2) {
		auto f = sp_compute(results.back(), net) * 20 - 10;
		results.push_back(sp_timestep(results.back(), f));
		results.push_back(sp_timestep(results.back(), 0));

		if (!is_upright(results.back())) break;
	}

	float f1 = 0.f;
	float f2 = 0.f;

	f1 = (float)t / N_T;

	if (10 * t > N_T) {
		float denom = 0;
		for (size_t i = 0, j = 0; i < results.size(); ++i) {
			if (is_upright(results[i]) && j++ > t - N_T / 10) {
				denom += fabsf((float)results[i].cart_x);
				denom += fabsf((float)results[i].cart_v);
				denom += fabsf((float)results[i].pole_θ);
				denom += fabsf((float)results[i].pole_v);
			}
		}

		f2 = 0.75f / denom;
	}

	return f1 * 0.1f + f2 * 0.9f;
}

// ===========================  DPV
// ===========================  DPV
// ===========================  DPV
// ===========================  DPV


struct Double_Pole_State {
	double cart_x = 0;
	double cart_v = 0;

	double pole1_θ = 0;
	double pole1_v = 0;
	double pole2_θ = 0;
	double pole2_v = 0;
};
auto dpv_compute = [](Double_Pole_State state, Network& net) -> float {
	float inputs[] = {
		(float)state.cart_x,
		(float)state.cart_v,
		(float)state.pole1_θ,
		(float)state.pole1_v,
		(float)state.pole2_θ,
		(float)state.pole2_v,
		1.f // biais
	};

	inputs[0] = inputs[0] / 2.4f;
	inputs[1] = sig(inputs[1]);
	inputs[2] = (atan2(sinf(inputs[2]), cosf(inputs[2]))) / (M_PI);
	inputs[3] = sig(inputs[3]);
	inputs[4] = (atan2(sinf(inputs[4]), cosf(inputs[4]))) / (M_PI);
	inputs[5] = sig(inputs[5]);

	float outputs[1];

	net.compute(inputs, 7, outputs, 1);

	return outputs[0];
};

auto is_upright_dp = [](Double_Pole_State x) {
	return
		std::abs(x.pole1_θ) < M_PI_4 &&
		std::abs(x.pole2_θ) < M_PI_4;
};

auto dp_compute = [](Double_Pole_State state, Network& net) -> float {
	float inputs[] = {
		(float)state.cart_x,
		(float)state.cart_v,
		(float)state.pole1_θ,
		(float)state.pole2_θ,
		1.f // biais
	};

	inputs[0] = inputs[0] / 2.4f;
	inputs[1] = sig(inputs[1]);
	inputs[2] = (atan2(sinf(inputs[2]), cosf(inputs[2]))) / (M_PI);
	inputs[3] = (atan2(sinf(inputs[4]), cosf(inputs[4]))) / (M_PI);

	float outputs[1];

	net.compute(inputs, 5, outputs, 1);

	return outputs[0];
};

auto dp_timestep = [](Double_Pole_State in, double f) -> Double_Pole_State {
	Double_Pole_State out = in;
	constexpr double m    = 0.1;     // mass of the pole
	constexpr double m2   = 0.01;    // mass of the pole2
	constexpr double m_c  = 1.0;     // mass of the cart
	constexpr double g    = 9.8;     // gravity
	constexpr double l    = 1.0;     // length of the pole
	constexpr double l2   = 0.1;     // length of the pole2
	constexpr double mu_c = 0.1;     // cart firction coeff
	constexpr double mu_p = 0.01;    // joint firction coeff
	constexpr double M    = m + m_c; // combined mass of the system
	constexpr double t_l  = 4.8;     // trach width in meter

	constexpr double dt = 0.01; // timestep in second

	// Extract state into named variables.
	// Precompute some reused values.
	double sin_theta = std::sin(in.pole1_θ);
	double cos_theta = std::cos(in.pole1_θ);
	double sin_theta2 = std::sin(in.pole2_θ);
	double cos_theta2 = std::cos(in.pole2_θ);
	double cos_theta_sqr = cos_theta * cos_theta;
	double thetav_sqr = in.pole1_v * in.pole1_v;
	double cos_theta2_sqr = cos_theta2 * cos_theta2;
	double thetav2_sqr = in.pole2_v * in.pole2_v;

	double l_hat = l / 2;
	double l2_hat = l2 / 2;

	// Calc cart horizontal acceleration.
	auto xa = (g * ((m*sin_theta*cos_theta) + (m2*sin_theta2*cos_theta2))
                - (7./3.) * (f +(m*l_hat*thetav_sqr*sin_theta) + (m2*l2_hat*thetav2_sqr*sin_theta2) - mu_c*in.cart_v)
                - (((mu_p*in.pole1_v*cos_theta)/l_hat) + ((mu_p*in.pole2_v*cos_theta2)/l2_hat)))
                / ((m*cos_theta_sqr) + (m2*cos_theta2_sqr) - (7./3.)*M);

	if (in.cart_x <= -t_l / 2 && xa < 0) xa = 0;
	if (in.cart_x >= +t_l / 2 && xa > 0) xa = 0;

	out.cart_x = in.cart_x + in.cart_v * dt;
	out.cart_v = in.cart_v + xa * dt;


	auto θa1 = (3./(7.*l_hat)) * (g*sin_theta - xa*cos_theta - ((mu_p * in.pole1_v)/(m*l_hat)));
	out.pole1_θ = in.pole1_θ + in.pole1_v * dt;
	out.pole1_v = in.pole1_v + θa1 * dt;


	auto θa2 = (3./(7.*l2_hat)) * (g*sin_theta2 - xa*cos_theta2 - ((mu_p * in.pole2_v)/(m2*l2_hat)));
	out.pole2_θ = in.pole2_θ + in.pole2_v * dt;
	out.pole2_v = in.pole2_v + θa2 * dt;

	if (out.cart_x <= -t_l / 2) {
		out.cart_x = -t_l / 2;
		out.cart_v = std::max(0.0, out.cart_v);
	}
	if (out.cart_x >= +t_l / 2) {
		out.cart_x = +t_l / 2;
		out.cart_v = std::min(0.0, out.cart_v);
	}

	return out;
};

void dpv_render(Network& net, void* user) noexcept {
	thread_local std::vector<Double_Pole_State> results;
	thread_local std::vector<double> actions;
	results.clear();
	actions.clear();

	Double_Pole_State state;
	state.cart_v = 0;
	state.cart_x = 0;
	state.pole1_v = 0;
	state.pole1_θ = (1.f / 180) * M_PI;
	state.pole2_v = 0;
	state.pole2_θ = 0;

	results.push_back(state);
	size_t N_T = 1000;
	for (size_t i = 0; i < N_T; i += 2) {
		auto f = dpv_compute(results.back(), net) * 20 - 10;
		results.push_back(dp_timestep(results.back(), f));
		results.push_back(dp_timestep(results.back(), 0));
		actions.push_back(f);
		actions.push_back(0);
	}

	thread_local float time = 0;
	thread_local bool run = false;

	if (ImGui::Button(run ? "Pause" : " Run ")) run = !run;
	ImGui::SameLine();
	ImGui::SliderFloat("Time", &time, 0, 1);

	if (run) {
		time += ImGui::GetIO().DeltaTime / 10;
	}

	time = fmodf(time, 1.f);
	time = std::max(time, 0.f);

	state = results[time * (results.size() - 1)];
	auto action = actions[time * (actions.size() - 1)];

	auto draw_list = ImGui::GetWindowDrawList();

	auto offset = ImGui::GetWindowPos();
	auto pos = ImGui::GetCursorPos();
	auto size = ImGui::GetContentRegionAvail();
	auto center = pos + size / 2;

	auto left_end = ImVec2{ pos.x + size.x / 20, pos.y + 4 * size.y / 5 };
	auto right_end = ImVec2{ pos.x + 19 * size.x / 20, pos.y + 4 * size.y / 5 };



	auto scale = (right_end.x - left_end.x) / 4.8f;

	auto velocity_height = 0.05;
	auto force_height = 0.035;

	auto t = ((float)state.cart_x + 2.4f) / 4.8f;
	auto cart_pos = ImVec2{ left_end.x + (right_end.x - left_end.x) * t, left_end.y };

	auto pole1_start = cart_pos;
	auto pole1_end = pole1_start + ImVec2{
		cosf(state.pole1_θ - M_PI_2) * 0.1f * scale,
		sinf(state.pole1_θ - M_PI_2) * 0.1f * scale
	};
	auto pole2_start = pole1_end;
	auto pole2_end = pole2_start + ImVec2{
		cosf(state.pole2_θ - M_PI_2) * 1.f * scale,
		sinf(state.pole2_θ - M_PI_2) * 1.f * scale
	};

	// draw left end
	draw_list->AddRectFilled(
		offset + left_end - ImVec2(0.05f * scale, 0.05f * scale),
		offset + left_end + ImVec2(0.05f * scale, 0.05f * scale),
		IM_COL32(255, 0, 0, 255)
	);
	// draw right end
	draw_list->AddRectFilled(
		offset + right_end - ImVec2(0.05f * scale, 0.05f * scale),
		offset + right_end + ImVec2(0.05f * scale, 0.05f * scale),
		IM_COL32(255, 0, 0, 255)
	);
	// draw base
	draw_list->AddRectFilled(
		offset + left_end  - ImVec2(0.0125f * scale, 0.0125f * scale),
		offset + right_end + ImVec2(0.0125f * scale, 0.0125f * scale),
		IM_COL32(255, 0, 0, 120)
	);
	// draw probe
	draw_list->AddRectFilled(
		offset + cart_pos + ImVec2(0.07f * scale, 0.07f * scale),
		offset + cart_pos - ImVec2(0.07f * scale, 0.07f * scale),
		IM_COL32(0, 255, 0, 255)
	);
	// draw poles
	draw_list->AddLine(
		offset + pole1_start,
		offset + pole1_end,
		IM_COL32(0, 0, 255, 200),
		0.02f * scale
	);
	draw_list->AddLine(
		offset + pole2_start,
		offset + pole2_end,
		IM_COL32(0, 0, 255, 200),
		0.02f * scale
	);
	// draw junction
	draw_list->AddCircle(offset + pole1_end, 0.05f * scale, IM_COL32(100, 100, 100, 100));

	draw_list->AddCircle(
		{ offset.x + cart_pos.x, offset.y + cart_pos.y },
		50,
		is_upright_dp(state) ? IM_COL32(0, 255, 0, 255) : IM_COL32(255, 0, 0, 255)
	);

	// Draw force.
	draw_list->AddRectFilled(
		{
			(float)(offset.x + cart_pos.x),
			(float)(offset.y + cart_pos.y + 0.25f * scale - force_height * scale)
		},
		{
			(float)(offset.x + cart_pos.x + action * scale),
			(float)(offset.y + cart_pos.y + 0.25f * scale + force_height * scale)
		},
		IM_COL32(200, 200, 200, 200)
	);

	// Draw velocity.
	draw_list->AddRectFilled(
		{
			(float)(offset.x + cart_pos.x),
			(float)(offset.y + cart_pos.y + 0.15f * scale - velocity_height * scale)
		},
		{
			(float)(offset.x + cart_pos.x + state.cart_v * scale),
			(float)(offset.y + cart_pos.y + 0.15f * scale + velocity_height * scale)
		},
		IM_COL32(100, 100, 100, 200)
	);

	ImGui::SetCursorPos(pos + size);
}

float dpv_fitness(Network& net, void*) noexcept {
	thread_local std::vector<Double_Pole_State> results;
	results.clear();

	Double_Pole_State state;
	state.cart_v = 0;
	state.cart_x = 0;
	state.pole1_v = 0;
	state.pole1_θ = (1.f / 180) * M_PI;
	state.pole2_v = 0;
	state.pole2_θ = 0;

	size_t N_T = 1000;

	results.push_back(state);

	size_t t = 0;
	for (; t < N_T; t += 2) {
		auto f = dpv_compute(results.back(), net) * 20 - 10;
		results.push_back(dp_timestep(results.back(), f));
		results.push_back(dp_timestep(results.back(), 0));

		if (!is_upright_dp(results.back())) break;
	}

	float f1 = 0.f;
	float f2 = 0.f;

	f1 = (float)t / N_T;

	if (10 * t > N_T) {
		float denom = 0;
		for (size_t i = 0, j = 0; i < results.size(); ++i) {
			if (is_upright_dp(results[i]) && j++ > t - N_T / 10) {
				denom += fabsf((float)results[i].cart_x);
				denom += fabsf((float)results[i].cart_v);
				denom += fabsf((float)results[i].pole1_θ);
				denom += fabsf((float)results[i].pole2_θ);
			}
		}

		f2 = 0.75f / denom;
	}

	return f1;
	return f1 * .1f + f2 * .9f;
}


void dp_render(Network& net, void* user) noexcept {
	thread_local std::vector<Double_Pole_State> results;
	thread_local std::vector<double> actions;
	results.clear();
	actions.clear();

	Double_Pole_State state;
	state.cart_v = 0;
	state.cart_x = 0;
	state.pole1_v = 0;
	state.pole1_θ = (1.f / 180) * M_PI;
	state.pole2_v = 0;
	state.pole2_θ = 0;

	results.push_back(state);
	auto is_upright = [](Double_Pole_State x) {
		return
			std::abs(x.pole1_θ) < M_PI_4 / 6 &&
			std::abs(x.pole2_θ) < M_PI_4 / 6;
	};

	size_t N_T = 1000;
	for (size_t i = 0; i < N_T; i += 2) {
		auto f = dp_compute(results.back(), net) * 20 - 10;
		results.push_back(dp_timestep(results.back(), f));
		results.push_back(dp_timestep(results.back(), 0));
		actions.push_back(f);
		actions.push_back(0);
	}

	thread_local float time = 0;
	thread_local bool run = false;
	if (run) {
		time += ImGui::GetIO().DeltaTime / 10;
	}

	if (ImGui::Button(run ? "Pause" : " Run ")) run = !run;
	ImGui::SameLine();
	ImGui::SliderFloat("Time", &time, 0, 1);

	time = std::min(1.f, std::max(time, 0.f));

	state = results[time * (results.size() - 1)];
	auto action = actions[time * (actions.size() - 1)];

	auto draw_list = ImGui::GetWindowDrawList();

	auto offset = ImGui::GetWindowPos();
	auto pos = ImGui::GetCursorPos();
	auto size = ImGui::GetContentRegionAvail();
	auto center = pos + size / 2;

	auto left_end = ImVec2{ pos.x + size.x / 20, pos.y + 4 * size.y / 5 };
	auto right_end = ImVec2{ pos.x + 19 * size.x / 20, pos.y + 4 * size.y / 5 };

	auto scale = (right_end.x - left_end.x) / 4.8f;

	auto velocity_height = 0.05;
	auto force_height = 0.035;

	auto t = ((float)state.cart_x + 2.4f) / 4.8f;
	auto cart_pos = ImVec2{ left_end.x + (right_end.x - left_end.x) * t, left_end.y };

	auto pole1_start = cart_pos;
	auto pole1_end = pole1_start + ImVec2{
		cosf(state.pole1_θ - M_PI_2) * 1.0f * scale,
		sinf(state.pole1_θ - M_PI_2) * 1.0f * scale
	};
	auto pole2_start = pole1_end;
	auto pole2_end = pole2_start + ImVec2{
		cosf(state.pole2_θ - M_PI_2) * 0.1f * scale,
		sinf(state.pole2_θ - M_PI_2) * 0.1f * scale
	};

	// draw left end
	draw_list->AddRectFilled(
		offset + left_end - ImVec2(0.05f * scale, 0.05f * scale),
		offset + left_end + ImVec2(0.05f * scale, 0.05f * scale),
		IM_COL32(255, 0, 0, 255)
	);
	// draw right end
	draw_list->AddRectFilled(
		offset + right_end - ImVec2(0.05f * scale, 0.05f * scale),
		offset + right_end + ImVec2(0.05f * scale, 0.05f * scale),
		IM_COL32(255, 0, 0, 255)
	);
	// draw base
	draw_list->AddRectFilled(
		offset + left_end  - ImVec2(0.0125f * scale, 0.0125f * scale),
		offset + right_end + ImVec2(0.0125f * scale, 0.0125f * scale),
		IM_COL32(255, 0, 0, 120)
	);
	// draw probe
	draw_list->AddRectFilled(
		offset + cart_pos + ImVec2(0.07f * scale, 0.07f * scale),
		offset + cart_pos - ImVec2(0.07f * scale, 0.07f * scale),
		IM_COL32(0, 255, 0, 255)
	);
	// draw poles
	draw_list->AddLine(
		offset + pole1_start,
		offset + pole1_end,
		IM_COL32(0, 0, 255, 200),
		0.02f * scale
	);
	draw_list->AddLine(
		offset + pole2_start,
		offset + pole2_end,
		IM_COL32(0, 0, 255, 200),
		0.02f * scale
	);
	// draw junction
	draw_list->AddCircle(offset + pole1_end, 0.05f * scale, IM_COL32(100, 100, 100, 100));

	draw_list->AddCircle(
		{ offset.x + cart_pos.x, offset.y + cart_pos.y },
		50,
		is_upright_dp(state) ? IM_COL32(0, 255, 0, 255) : IM_COL32(255, 0, 0, 255)
	);

	// Draw force.
	draw_list->AddRectFilled(
		{
			(float)(offset.x + cart_pos.x),
			(float)(offset.y + cart_pos.y + 0.25f * scale - force_height * scale)
		},
		{
			(float)(offset.x + cart_pos.x + action * scale),
			(float)(offset.y + cart_pos.y + 0.25f * scale + force_height * scale)
		},
		IM_COL32(200, 200, 200, 200)
	);

	// Draw velocity.
	draw_list->AddRectFilled(
		{
			(float)(offset.x + cart_pos.x),
			(float)(offset.y + cart_pos.y + 0.15f * scale - velocity_height * scale)
		},
		{
			(float)(offset.x + cart_pos.x + state.cart_v * scale),
			(float)(offset.y + cart_pos.y + 0.15f * scale + velocity_height * scale)
		},
		IM_COL32(100, 100, 100, 200)
	);

	ImGui::SetCursorPos(pos + size);
}

float dp_fitness(Network& net, void*) noexcept {
	thread_local std::vector<Double_Pole_State> results;
	results.clear();

	Double_Pole_State state;
	state.cart_v = 0;
	state.cart_x = 0;
	state.pole1_v = 0;
	state.pole1_θ = (1.f / 180) * M_PI;
	state.pole2_v = 0;
	state.pole2_θ = 0;

	size_t N_T = 1000;

	results.push_back(state);
	size_t t = 0;
	for (; t < N_T; t += 2) {
		auto f = dp_compute(results.back(), net) * 20 - 10;
		results.push_back(dp_timestep(results.back(), f));
		results.push_back(dp_timestep(results.back(), 0));

		if (!is_upright_dp(results.back())) break;
	}

	float f1 = 0.f;
	float f2 = 0.f;

	f1 = (float)t / N_T;

	if (10 * t > N_T) {
		float denom = 0;
		for (size_t i = 0, j = 0; i < results.size(); ++i) {
			if (is_upright_dp(results[i]) && j++ > t - N_T / 10) {
				denom += fabsf((float)results[i].cart_x);
				denom += fabsf((float)results[i].cart_v);
				denom += fabsf((float)results[i].pole1_θ);
				denom += fabsf((float)results[i].pole2_θ);
			}
		}

		f2 = 0.75f / denom;
	}

	return f1 * 0.1f + f2 * 0.9f;
}

void hdpv_render(Network& net, void*) noexcept {
	// compute CPPN.
	Substrate input;
	Substrate hidden;
	Substrate output;
	output.density = 1;

	mark_node_plan(input , 1, 1, Network::Node_Kind::Sensor); // amount of negative cart velocity
	mark_node_plan(input , 2, 1, Network::Node_Kind::Sensor); // cart position
	mark_node_plan(input , 3, 1, Network::Node_Kind::Sensor); // amount of positive cart velocity
	mark_node_plan(input , 1, 2, Network::Node_Kind::Sensor); // amount of negative pole1 velocity
	mark_node_plan(input , 2, 2, Network::Node_Kind::Sensor); // pole1 position
	mark_node_plan(input , 3, 2, Network::Node_Kind::Sensor); // amount of positive pole1 velocity
	mark_node_plan(input , 1, 3, Network::Node_Kind::Sensor); // amount of negative pole2 velocity
	mark_node_plan(input , 2, 3, Network::Node_Kind::Sensor); // pole2 position
	mark_node_plan(input , 3, 3, Network::Node_Kind::Sensor); // amount of positive pole2 velocity
	mark_node_plan(output, 0, 0, Network::Node_Kind::Output);

	Network hyper_net = hyper({input, hidden, output}, net);

	return dpv_render(net, nullptr);
}

float hdpv_fitness(Network& net, void*) noexcept {
	// compute CPPN.
	Substrate input;
	Substrate hidden;
	Substrate output;
	output.density = 1;

	mark_node_plan(input , 1, 1, Network::Node_Kind::Sensor); // amount of negative cart velocity
	mark_node_plan(input , 2, 1, Network::Node_Kind::Sensor); // cart position
	mark_node_plan(input , 3, 1, Network::Node_Kind::Sensor); // amount of positive cart velocity
	mark_node_plan(input , 1, 2, Network::Node_Kind::Sensor); // amount of negative pole1 velocity
	mark_node_plan(input , 2, 2, Network::Node_Kind::Sensor); // pole1 position
	mark_node_plan(input , 3, 2, Network::Node_Kind::Sensor); // amount of positive pole1 velocity
	mark_node_plan(input , 1, 3, Network::Node_Kind::Sensor); // amount of negative pole2 velocity
	mark_node_plan(input , 2, 3, Network::Node_Kind::Sensor); // pole2 position
	mark_node_plan(input , 3, 3, Network::Node_Kind::Sensor); // amount of positive pole2 velocity
	mark_node_plan(output, 0, 0, Network::Node_Kind::Output);


	Network hyper_net = hyper({input, hidden, output}, net);

	return dpv_fitness(hyper_net, nullptr);
}