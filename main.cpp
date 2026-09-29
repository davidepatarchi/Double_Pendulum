#include "pendulum.hpp"
#include "render.hpp"

#include <SFML/Graphics.hpp>
#include <iostream>
#include <numbers>

int main() {
  try {
    sf::RenderWindow window(sf::VideoMode({800, 800}), "Double Pendulum");

    float lenght1 = 150.f;
    double mass1 = 10.;
    float theta1 = 1.5f;

    float lenght2 = 150.f;
    double mass2 = 10.;
    float theta2 = 0.;

    sf::Vector2f origin(400.f, 400.f);

    pnd::Pendulum pendulum1{lenght1 / pix::scale, mass1, theta1};
    pnd::Pendulum pendulum2{lenght2 / pix::scale, mass2, theta2};

    while (window.isOpen()) {
      while (const auto event = window.pollEvent()) {
        if (event->is<sf::Event::Closed>())
          window.close();
      }

      pendulum1.evolution(pendulum2);

      float x1 = float(pendulum1.pos_x(0.));
      float x2 = float(pendulum2.pos_x(x1));
      float y1 = float(pendulum1.pos_y(0.));
      float y2 = float(pendulum2.pos_y(x2));

      sf::RectangleShape rod1;
      rod1.setSize({lenght1, 4.f});
      rod1.setFillColor(sf::Color::White);
      rod1.setOrigin({0.f, 2.f});
      rod1.setPosition(origin);
      rod1.setRotation(
          sf::radians(float((std::numbers::pi / 2.0) - pendulum1.theta())));

      sf::RectangleShape rod2;
      rod1.setSize({lenght2, 4.f});
      rod1.setFillColor(sf::Color::White);
      rod1.setOrigin({0.f, 2.f});
      rod1.setPosition(sf::Vector2f({x1, y1}));
      rod1.setRotation(
          sf::radians(float((std::numbers::pi / 2.0) - pendulum2.theta())));

      sf::CircleShape point1(10.f);
      point1.setFillColor(sf::Color::Red);
      point1.setOrigin({10.f, 10.f});
      point1.setPosition({pix_pos_x(x1), pix_pos_y(y1)});

      sf::CircleShape point2(10.f);
      point1.setFillColor(sf::Color::Red);
      point1.setOrigin({10.f, 10.f});
      point1.setPosition({pix_pos_x(x2), pix_pos_y(y2)});

      sf::CircleShape pivot(8.f);
      pivot.setFillColor(sf::Color::White);
      pivot.setOrigin({8.f, 8.f});
      pivot.setPosition(origin);

      window.clear(sf::Color::Black);

      window.draw(rod1);
      window.draw(point1);
      window.draw(rod2);
      window.draw(point2);
      window.draw(pivot);

      window.display();
    }
  } catch (const std::invalid_argument &error) {
    std::cerr << error.what();
    return EXIT_FAILURE;
  }
}
