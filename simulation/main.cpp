#include "../rendering/render.hpp"
#include "pendulum.hpp"

#include <fstream>
#include <iostream>

int main(int argc, char *argv[]) {
  if (argc > 3) {
    throw std::runtime_error("Error: expected not more than 3 imput values");
  }

  sf::RenderWindow window(sf::VideoMode({800, 800}), "Double Pendulum");
  window.setFramerateLimit(60);

  float lenght1 = 150.f;
  double mass1 = 10.;
  float lenght2 = 150.f;
  double mass2 = 10.;

  double theta1;
  if (argc > 1) {
    theta1 = std::stod(argv[1]) * std::numbers::pi / 180.;
  } else {
    theta1 = 0.f;
  }
  double theta2;
  if (argc > 2) {
    theta2 = std::stod(argv[2]) * std::numbers::pi / 180.;
  } else {
    theta2 = 0.f;
  }

  pnd::State state{theta1, theta2, 0., 0.};
  pnd::DoublePendulum pend{lenght1 / pix::scale, mass1, lenght2 / pix::scale,
                           mass2, state};

  sf::Vector2f origin(400.f, 400.f);

  sf::RectangleShape vert1;
  vert1.setSize({2.f, lenght1 / 3});
  vert1.setFillColor(sf::Color::Magenta);
  vert1.setOrigin({0., 1.f});
  vert1.setPosition(origin);

  bool dragging1 = false;
  bool dragging2 = false;

  double T{0.};
  double K{0.};
  double U{0.};
  double E{0.};

  std::ofstream file("../data/energy.csv");
  if (!file.is_open()) {
    std::cerr << "Could not open energy.csv\n";
  }
  file << "Time,Kinetic,Potential,Total\n";

  while (window.isOpen()) {
    while (const auto event = window.pollEvent()) {

      if (event->is<sf::Event::Closed>()) {
        window.close();
      }

      if (const auto *key = event->getIf<sf::Event::KeyPressed>()) {
        if (key->code == sf::Keyboard::Key::Space) {
          pend.state() = {0., 0., 0., 0.};
        }
      }

      if (const auto *mouse = event->getIf<sf::Event::MouseButtonPressed>()) {
        if (mouse->button == sf::Mouse::Button::Left) {

          sf::Vector2f mouse_position{float(mouse->position.x),
                                      float(mouse->position.y)};

          sf::Vector2f position1{pix_pos_x(pend.pos_x1(0.)),
                                 pix_pos_y(pend.pos_y1(0.))};

          sf::Vector2f position2{pix_pos_x(pend.pos_x2(pend.pos_x1(0.))),
                                 pix_pos_y(pend.pos_y2(pend.pos_y1(0.)))};

          float distance1 = std::hypot(mouse_position.x - position1.x,
                                       mouse_position.y - position1.y);

          float distance2 = std::hypot(mouse_position.x - position2.x,
                                       mouse_position.y - position2.y);

          if (distance1 < 20.f) {
            dragging1 = true;
          } else if (distance2 < 20.f) {
            dragging2 = true;
          }
        }
      }

      if (const auto *mouse = event->getIf<sf::Event::MouseButtonReleased>()) {
        if (mouse->button == sf::Mouse::Button::Left) {
          dragging1 = false;
          dragging2 = false;
        }
      }
    }

    if (dragging1 || dragging2) {
      sf::Vector2i mouse_pixel = sf::Mouse::getPosition(window);
      sf::Vector2f mouse{float(mouse_pixel.x), float(mouse_pixel.y)};

      if (dragging1) {
        float dx = mouse.x - origin.x;
        float dy = mouse.y - origin.y;

        pend.state().theta1 = std::atan2(dx, dy);
        pend.state().omega1 = 0.0;
      }

      if (dragging2) {

        float x1 = float(pend.pos_x1(0.));
        float y1 = float(pend.pos_y1(0.));

        sf::Vector2f pivot{pix_pos_x(x1), pix_pos_y(y1)};

        float dx = mouse.x - pivot.x;
        float dy = mouse.y - pivot.y;

        pend.state().theta2 = std::atan2(dx, dy);
        pend.state().omega2 = 0.0;
      }
    }

    for (int i{0}; i < 17; ++i) {
      pend.evolution();
    }

    float x1 = float(pend.pos_x1(0.));
    float x2 = float(pend.pos_x2(x1));
    float y1 = float(pend.pos_y1(0.));
    float y2 = float(pend.pos_y2(y1));

    sf::RectangleShape rod1;
    rod1.setSize({4.f, lenght1});
    rod1.setFillColor(sf::Color::Magenta);
    rod1.setOrigin({2.f, 0.f});
    rod1.setPosition(origin);
    rod1.setRotation(sf::radians(float(-pend.state().theta1)));

    sf::RectangleShape rod2;
    rod2.setSize({4.f, lenght2});
    rod2.setFillColor(sf::Color::Cyan);
    rod2.setOrigin({2.f, 0.f});
    rod2.setPosition(sf::Vector2f({pix_pos_x(x1), pix_pos_y(y1)}));
    rod2.setRotation(sf::radians(float(-pend.state().theta2)));

    sf::CircleShape point1(10.f);
    point1.setFillColor(sf::Color::White);
    point1.setOrigin({10.f, 10.f});
    point1.setPosition({pix_pos_x(x1), pix_pos_y(y1)});

    sf::CircleShape point2(10.f);
    point2.setFillColor(sf::Color::White);
    point2.setOrigin({10.f, 10.f});
    point2.setPosition({pix_pos_x(x2), pix_pos_y(y2)});

    sf::CircleShape pivot(5.f);
    pivot.setFillColor(sf::Color::White);
    pivot.setOrigin({5.f, 5.f});
    pivot.setPosition(origin);

    sf::RectangleShape vert2;
    vert2.setSize({2.f, lenght2 / 3});
    vert2.setFillColor(sf::Color::Cyan);
    vert2.setOrigin({1.f, 0.f});
    vert2.setPosition(sf::Vector2f{pix_pos_x(x1), pix_pos_y(y1)});

    sf::VertexArray arc1 =
        drawArc(origin, lenght1 / 4.f, float(std::numbers::pi / 2.f),
                float(std::numbers::pi / 2.f - pend.state().theta1),
                sf::Color::Magenta);

    sf::VertexArray arc2 = drawArc(
        sf::Vector2f{pix_pos_x(x1), pix_pos_y(y1)}, lenght1 / 4.f,
        float(std::numbers::pi / 2.f),
        float(std::numbers::pi / 2.f - pend.state().theta2), sf::Color::Cyan);

    T += 17 * pnd::constants::dt;
    K = pend.kinetic();
    U = pend.potential();
    E = K + U;
    file << T << ',' << K << ',' << U << ',' << E << '\n';

    window.clear(sf::Color::Black);

    window.draw(rod1);
    window.draw(rod2);
    window.draw(vert1);
    window.draw(vert2);
    window.draw(arc1);
    window.draw(arc2);
    window.draw(point1);
    window.draw(point2);
    window.draw(pivot);

    window.display();
  }
  file.close();
}