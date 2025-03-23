#pragma once
#include <algorithm>
#include <ranges>
#include <tuple>

#include "Expr.h"
#include "ExprFactory.h"

namespace Expression {

inline const auto unity = ExprFactory::number(1.0);
inline const auto zero = ExprFactory::number(0.0);

// Template to detect std::variant
template <typename T>
concept VariantType = requires {
  typename std::variant_size<std::remove_cvref_t<T>>::type;
  requires std::variant_size_v<std::remove_cvref_t<T>> > 0;
};

// Helper templates for overloaded visitors
template <class... Ts>
struct overloaded : Ts... {
  using Ts::operator()...;
};
template <class... Ts>
overloaded(Ts...) -> overloaded<Ts...>;

// Helper that captures variants
template <VariantType... Variants>
class MatchMaker {
 public:
  explicit MatchMaker(const Variants&... vars) : variants_(vars...) {}

  // With function to capture the lambdas
  template <typename... Lambdas>
  auto with(Lambdas&&... lambdas) const {
    return std::apply(
        [&](const auto&... vars) {
          return std::visit(overloaded{std::forward<Lambdas>(lambdas)...},
                            vars...);
        },
        variants_);
  }

 private:
  std::tuple<const Variants&...> variants_;
};

template <VariantType... Variants>
auto match(const Variants&... variants) {
  return MatchMaker<Variants...>(variants...);
}

inline auto match(const Expr& expr) { return match(expr.getImpl()); }

// Helper template for static_assert that depends on a type
template <typename T>
inline constexpr bool always_false_v = false;

// Helper templates for type checking
template <typename T>
concept NamedNullaryType = std::is_base_of_v<NamedNullaryExpr, std::decay_t<T>>;

template <typename T>
concept UnaryType = std::is_base_of_v<UnaryExpr, std::decay_t<T>>;

template <typename T>
concept NaryType = std::is_base_of_v<NaryExpr, std::decay_t<T>>;

template <typename T, VariantType Variant>
bool is(const Variant& v) {
  return std::holds_alternative<T>(v);
}

template <typename T>
bool is(const Expr& expr) {
  return is<T>(expr.getImpl());
}

template <typename TLambda>
std::vector<Expr> transform(const std::vector<Expr>& terms,
                            const TLambda& lambda) {
  std::vector<Expr> transformedTerms;
  transformedTerms.reserve(terms.size());
  std::ranges::transform(terms, std::back_inserter(transformedTerms),
                         [&lambda](const auto& t) { return lambda(t); });
  return transformedTerms;
}

}  // namespace Expression