#include "pendulum.hpp"
#include "render.hpp"

#include <SFML/Graphics.hpp>
#include <iostream>
#include <numbers>

int main(int argc, char *argv[]) {
  if (argc != 3) {
    throw std::runtime_error("Error: expected 3 imput values");
  }

  sf::RenderWindow window(sf::VideoMode({800, 800}), "Double Pendulum");

  float lenght1 = 150.f;
  double mass1 = 10.;
  float theta1 = std::stof(argv[1]) * 0.017f;

  float lenght2 = 150.f;
  double mass2 = 10.;
  float theta2 = std::stof(argv[2]) * 0.017f;

  sf::Vector2f origin(400.f, 400.f);

  pnd::Pendulum pendulum1{lenght1 / pix::scale, mass1, theta1};
  pnd::Pendulum pendulum2{lenght2 / pix::scale, mass2, theta2};

  sf::RectangleShape vert1;
  vert1.setSize({2.f, lenght1 / 3});
  vert1.setFillColor(sf::Color::Yellow);
  vert1.setOrigin({0., 1.f});
  vert1.setPosition(origin);

  while (window.isOpen()) {
    while (const auto event = window.pollEvent()) {
      if (event->is<sf::Event::Closed>())
        window.close();
    }

    pendulum1.evolution(pendulum2);

    float x1 = float(pendulum1.pos_x(0.));
    float x2 = float(pendulum2.pos_x(x1));
    float y1 = float(pendulum1.pos_y(0.));
    float y2 = float(pendulum2.pos_y(y1));

    sf::RectangleShape rod1;
    rod1.setSize({4.f, lenght1});
    rod1.setFillColor(sf::Color::White);
    rod1.setOrigin({2.f, 0.f});
    rod1.setPosition(origin);
    rod1.setRotation(sf::radians(float(-pendulum1.theta())));

    sf::RectangleShape rod2;
    rod2.setSize({4.f, lenght2});
    rod2.setFillColor(sf::Color::White);
    rod2.setOrigin({2.f, 0.f});
    rod2.setPosition(sf::Vector2f({pix_pos_x(x1), pix_pos_y(y1)}));
    rod2.setRotation(sf::radians(float(-pendulum2.theta())));

    sf::CircleShape point1(10.f);
    point1.setFillColor(sf::Color::Red);
    point1.setOrigin({10.f, 10.f});
    point1.setPosition({pix_pos_x(x1), pix_pos_y(y1)});

    sf::CircleShape point2(10.f);
    point2.setFillColor(sf::Color::Red);
    point2.setOrigin({10.f, 10.f});
    point2.setPosition({pix_pos_x(x2), pix_pos_y(y2)});

    sf::CircleShape pivot(8.f);
    pivot.setFillColor(sf::Color::White);
    pivot.setOrigin({8.f, 8.f});
    pivot.setPosition(origin);

    sf::RectangleShape vert2;
    vert2.setSize({2.f, lenght2 / 3});
    vert2.setFillColor(sf::Color::Yellow);
    vert2.setOrigin({1.f, 0.f});
    vert2.setPosition(sf::Vector2f{pix_pos_x(x1), pix_pos_y(y1)});

    window.clear(sf::Color::Black);

    window.draw(rod1);
    window.draw(point1);
    window.draw(rod2);
    window.draw(point2);
    window.draw(pivot);
    window.draw(vert1);
    window.draw(vert2);

    window.display();
  }
}
