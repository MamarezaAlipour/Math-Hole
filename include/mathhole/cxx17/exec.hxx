// Copyright (c) 2022 Geek pww
// All rights reserved

#ifndef MATHHOLE_CXX17_EXEC_HXX
#define MATHHOLE_CXX17_EXEC_HXX

#include <type_traits>

namespace mathhole {
	namespace exec {

		//----------------------------------------------------------------------------//
		// Execution policy classes

		struct sequential_policy {};
		struct parallel_policy {};
		struct parallel_unsequential_policy {};
		struct unsequential_policy {};

		//----------------------------------------------------------------------------//
		// Execution policy instances

		constexpr auto seq = sequential_policy {};
		constexpr auto par = parallel_policy {};
		constexpr auto par_unseq = parallel_unsequential_policy {};
		constexpr auto unseq = unsequential_policy {};

		//----------------------------------------------------------------------------//
		// Execution policy type support

		// Type category

		template <typename T>
		struct is_execution_policy : std::false_type {};

		template <typename T>
		inline constexpr bool is_execution_policy_v = is_execution_policy<T>::value;

		// Type category specializations

		template <>
		struct is_execution_policy<sequential_policy> : std::true_type {};

		template <>
		struct is_execution_policy<parallel_policy> : std::true_type {};

		template <>
		struct is_execution_policy<parallel_unsequential_policy> : std::true_type {};

		template <>
		struct is_execution_policy<unsequential_policy> : std::true_type {};

	} // namespace exec
} // namespace mathhole

#endif