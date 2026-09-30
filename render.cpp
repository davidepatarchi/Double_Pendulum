#include "render.hpp"

float pix_pos_x(double x) { return float(400 + x * pix::scale); }
float pix_pos_y(double y) { return float(400 - y * pix::scale); }

sf::VertexArray drawArc(sf::Vector2f center, float radius, float start_angle,
                        float end_angle, sf::Color color) {

  sf::VertexArray arc(sf::PrimitiveType::LineStrip);

  constexpr int points = 100;

  for (int i = 0; i <= points; ++i) {

    float t = static_cast<float>(i) / points;
    float angle = start_angle + t * (end_angle - start_angle);

    sf::Vector2f point{center.x + radius * std::cos(angle),
                       center.y + radius * std::sin(angle)};

    arc.append(sf::Vertex{point, color});
  }

  return arc;
}