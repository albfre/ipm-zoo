#pragma once
#include <algorithm>
#include <ranges>
#include <tuple>

#include "Expression.h"

namespace Expression {

inline const auto unity = ExprFactory::number(1.0);
inline const auto zero = ExprFactory::number(0.0);

// Template to detect std::variant
template <typename T>
struct is_variant : std::false_type {};

template <typename... Types>
struct is_variant<std::variant<Types...>> : std::true_type {};

template <typename T>
inline constexpr bool is_variant_v = is_variant<T>::value;

template <typename T>
concept VariantType =
    requires { requires is_variant_v<std::remove_cvref_t<T>>; };

// Helper templates for overloaded visitors
template <class... Ts>
struct overloaded : Ts... {
  using Ts::operator()...;
};
template <class... Ts>
overloaded(Ts...) -> overloaded<Ts...>;

template <VariantType Variant, typename... Lambdas>
auto match(const Variant& var, Lambdas&&... lambdas) {
  return std::visit(overloaded{std::forward<Lambdas>(lambdas)...}, var);
}

template <typename... Lambdas>
auto match(const Expr& var, Lambdas&&... lambdas) {
  return match(var.getImpl(), std::forward<Lambdas>(lambdas)...);
}

// Helper template for static_assert that depends on a type
template <typename T>
inline constexpr bool always_false_v = false;

// Helper templates for type checking
template <typename T>
struct is_named_nullary : std::is_base_of<NamedNullaryExpr, std::decay_t<T>> {};

template <typename T>
inline constexpr bool is_named_nullary_v = is_named_nullary<T>::value;

template <typename T>
struct is_unary : std::is_base_of<UnaryExpr, std::decay_t<T>> {};

template <typename T>
inline constexpr bool is_unary_v = is_unary<T>::value;

template <typename T>
struct is_nary : std::is_base_of<NaryExpr, std::decay_t<T>> {};

template <typename T>
inline constexpr bool is_nary_v = is_nary<T>::value;

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