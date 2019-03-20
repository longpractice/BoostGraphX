#pragma once

namespace bglx {

	//Methods like boost::vertices, boost::in_edges return pair
	//of begin iterator and end iterator. It is not syntax sweet for ranged
	//for loops, we use this wrapper to facilitate ranged for loops
	//these iterators are normally copiable
	template <typename TIterator>
	struct Range {
		Range(const std::pair<TIterator, TIterator>& beginEnd)
			: m_begin{ beginEnd.first }
		, m_end{ beginEnd.second } {};

		TIterator m_begin;
		TIterator m_end;
		TIterator begin() const
		{
			return m_begin;
		}
		TIterator end() const
		{
			return m_end;
		}
		bool empty() const
		{
			return m_begin == m_end;
		}
	};
}
