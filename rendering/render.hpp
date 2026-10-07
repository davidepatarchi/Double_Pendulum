#include <SFML/Graphics.hpp>

namespace pix {
inline constexpr int scale = 300; // 1 m is 400 pixel
}

float pix_pos_x(double x);
float pix_pos_y(double y);

sf::VertexArray drawArc(sf::Vector2f center, float radius, float start_angle,
                        float end_angle, sf::Color color);