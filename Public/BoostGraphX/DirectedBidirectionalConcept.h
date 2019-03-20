#pragma once
#include <boost/graph/graph_concepts.hpp>

namespace bglx
{
	//
	// In many cases, we require the graph to be a directed and bidirectional graph
	// The bidirectional requirement will speed up most of the algorithms
	// It however, could add some space cost

	// A future improvement would be to generically support non bidirectional graph
	template <typename G>
	struct DirectedBidirectionalGraphConcept {
		BOOST_CONCEPT_ASSERT((boost::BidirectionalGraphConcept<G>));
		static_assert(std::is_convertible<
			boost::graph_traits<G>::directed_category, boost::directed_tag>::value);
	};
}