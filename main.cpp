#include "pendulum.hpp"
#include "render.hpp"

#include <SFML/Graphics.hpp>
#include <iostream>

int main() {
  try {
    sf::RenderWindow window(sf::VideoMode({800, 800}), "Double Pendulum");

    float length = 300.f;
    float angle = 0.5f;
    sf::Vector2f origin(400.f, 200.f);
    double m = 1.;

    pnd::Pendulum pend({origin.x / pix::scale, origin.y / pix::scale},
                       length / pix::scale, angle, m);

    while (window.isOpen()) {
      while (const auto event = window.pollEvent()) {
        if (event->is<sf::Event::Closed>())
          window.close();
      }

      pend.evolution();

      sf::RectangleShape rod;
      rod.setSize({length, 4.f});
      rod.setFillColor(sf::Color::White);
      rod.setOrigin({0.f, 2.f});
      rod.setPosition(origin);
      rod.setRotation(sf::radians(angle));

      sf::CircleShape mass(20.f);
      mass.setFillColor(sf::Color::Red);
      mass.setOrigin({20.f, 20.f});
      mass.setPosition({float(pend.position().x * pix::scale),
                        float(pend.position().y * pix::scale)});

      sf::CircleShape pivot(8.f);
      pivot.setFillColor(sf::Color::White);
      pivot.setOrigin({8.f, 8.f});
      pivot.setPosition(origin);

      window.clear(sf::Color::Black);

      window.draw(rod);
      window.draw(mass);
      window.draw(pivot);

      window.display();
    }
  } catch (const std::invalid_argument &error) {
    std::cerr << error.what();
    return EXIT_FAILURE;
  }
}
