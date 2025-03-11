#pragma once
#include <tuple>
#include "Expression.h"

namespace Expression {

inline const auto unity = ExprFactory::number(1.0);
inline const auto zero = ExprFactory::number(0.0);

// Helper templates for overloaded visitors
template <class... Ts>
struct overloaded : Ts... {
  using Ts::operator()...;
};
template <class... Ts>
overloaded(Ts...) -> overloaded<Ts...>;

template <typename Variant, typename... Lambdas>
auto match(Variant&& var, Lambdas&&... lambdas) {
  return std::visit(overloaded{std::forward<Lambdas>(lambdas)...},
                    std::forward<Variant>(var));
}

template <typename Variant, typename... Lambdas>
void matchDefaultDoNothing(Variant&& var, Lambdas&&... lambdas) {
  std::visit(overloaded{[](const auto&) {}, std::forward<Lambdas>(lambdas)...},
             std::forward<Variant>(var));
}

// Helper template for static_assert that depends on a type
template <typename T>
inline constexpr bool always_false_v = false;

// Helper templates for type checking
template <typename T, typename Tuple>
struct is_any_of;

template <typename T, typename... Types>
struct is_any_of<T, std::tuple<Types...>>
    : std::disjunction<std::is_same<T, Types>...> {};

template <typename T, typename Tuple>
inline constexpr bool is_any_of_v = is_any_of<T, Tuple>::value;

template <typename T>
bool is(const Expr& expr) {
  return std::holds_alternative<T>(expr.getImpl());
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