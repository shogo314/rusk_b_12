#include <numeric>

using RuskSuffixIndices = std::vector<std::size_t>;

inline Obj<RuskVec<std::int64_t>> rusk_suffix_indices(const RuskSuffixIndices& indices) {
    auto result = make_obj<RuskVec<std::int64_t>>(RuskVec<std::int64_t>{});
    result->reserve(indices.size());
    for (auto index : indices) {
        result->push_back(make_obj<std::int64_t>(static_cast<std::int64_t>(index)));
    }
    return result;
}

template<class Read>
Obj<RuskVec<std::int64_t>> rusk_suffix_array(std::size_t n, Read read) {
    RuskSuffixIndices order(n);
    std::iota(order.begin(), order.end(), std::size_t{0});
    if (n == 0) return rusk_suffix_indices(order);
    std::sort(order.begin(), order.end(), [&](auto a, auto b) { return read(a) < read(b); });

    RuskSuffixIndices rank(n), next_rank(n), by_tail(n), offsets(n);
    std::size_t classes = 1;
    for (std::size_t i = 1; i < n; ++i) {
        if (read(order[i - 1]) != read(order[i])) ++classes;
        rank[order[i]] = classes - 1;
    }

    for (std::size_t width = 1; classes < n; width *= 2) {
        // Missing tails precede all existing tails. The rest follow the previous suffix order.
        std::size_t output = 0;
        for (std::size_t i = n - width; i < n; ++i) by_tail[output++] = i;
        for (auto i : order) {
            if (i >= width) by_tail[output++] = i - width;
        }

        std::fill_n(offsets.begin(), classes, std::size_t{0});
        for (auto i : by_tail) ++offsets[rank[i]];
        std::size_t position = 0;
        for (std::size_t i = 0; i < classes; ++i) {
            const auto count = offsets[i];
            offsets[i] = position;
            position += count;
        }
        for (auto i : by_tail) order[offsets[rank[i]]++] = i;

        const auto tail_rank = [&](auto i) { return i < n - width ? rank[i + width] + 1 : 0; };
        classes = 1;
        next_rank[order[0]] = 0;
        for (std::size_t i = 1; i < n; ++i) {
            const auto a = order[i - 1], b = order[i];
            if (rank[a] != rank[b] || tail_rank(a) != tail_rank(b)) ++classes;
            next_rank[b] = classes - 1;
        }
        rank.swap(next_rank);
    }
    return rusk_suffix_indices(order);
}

template<class Read>
Obj<RuskVec<std::int64_t>> rusk_lcp_array(
        std::size_t n, Read read, const RuskVec<std::int64_t>& suffixes, Loc loc) {
    if (n != suffixes.size()) fail(loc, "invalid suffix array length");
    RuskSuffixIndices inverse(n, n), order(n), lcp(n == 0 ? 0 : n - 1);
    for (std::size_t i = 0; i < n; ++i) {
        const auto start = *suffixes[i];
        if (start < 0 || static_cast<std::uint64_t>(start) >= n || inverse[start] != n) {
            fail(loc, "invalid suffix array index");
        }
        order[i] = static_cast<std::size_t>(start);
        inverse[order[i]] = i;
    }
    const auto tail_rank = [&](auto i) { return i + 1 < n ? inverse[i + 1] + 1 : 0; };
    for (std::size_t i = 1; i < n; ++i) {
        const auto a = order[i - 1], b = order[i];
        if (read(a) > read(b) || (read(a) == read(b) && tail_rank(a) >= tail_rank(b))) {
            fail(loc, "suffix array is not sorted");
        }
    }

    std::size_t common = 0;
    for (std::size_t i = 0; i < n; ++i) {
        if (inverse[i] == 0) {
            common = 0;
            continue;
        }
        const auto previous = order[inverse[i] - 1];
        while (common < n - i && common < n - previous && read(i + common) == read(previous + common)) {
            ++common;
        }
        lcp[inverse[i] - 1] = common;
        if (common != 0) --common;
    }
    return rusk_suffix_indices(lcp);
}
