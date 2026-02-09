#include <SFML/Graphics.hpp>
#include <cmath>
#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <vector> // Добавил для удобного хранения хвоста
#include <optional>

// --- КОНСТАНТЫ ---
const float G = 9.81f;
const float M1 = 1.0f;
const float M2 = 1.0f;
const float L1 = 1.0f;
const float L2 = 1.0f;
const float SCALE = 100.0f;
const float OFFSET_X = 400.0f;
const float OFFSET_Y = 200.0f;

// Структура состояния
struct State {
    float a1 = 0; float a2 = 0;
    float w1 = 0; float w2 = 0;
};

// --- 1. ДИФФЕРЕНЦИАЛЬНЫЕ УРАВНЕНИЯ ---
State getDerivatives(const State& s) {
    State d;
    d.a1 = s.w1;
    d.a2 = s.w2;

    float del = s.a1 - s.a2;
    float den1 = (M1 + M2) * L1 - M2 * L1 * cos(del) * cos(del);
    float den2 = (L2 / L1) * den1;

    float num1 = -G * (2 * M1 + M2) * sin(s.a1)
        - M2 * G * sin(s.a1 - 2 * s.a2)
        - 2 * sin(del) * M2 * (s.w2 * s.w2 * L2 + s.w1 * s.w1 * L1 * cos(del));

    float num2 = 2 * sin(del) * (s.w1 * s.w1 * L1 * (M1 + M2) + G * (M1 + M2) * cos(s.a1)
        + s.w2 * s.w2 * L2 * M2 * cos(del));

    d.w1 = num1 / den1;
    d.w2 = num2 / den2;
    return d;
}

// --- 2. ЧИСЛЕННЫЕ МЕТОДЫ ---
State solveEuler(State s, float dt) {
    State d = getDerivatives(s);
    s.a1 += d.a1 * dt; s.a2 += d.a2 * dt;
    s.w1 += d.w1 * dt; s.w2 += d.w2 * dt;
    return s;
}

State solveTrapezoidal(State s, float dt) {
    State d1 = getDerivatives(s);
    State s_pred = s;
    s_pred.a1 += d1.a1 * dt; s_pred.a2 += d1.a2 * dt;
    s_pred.w1 += d1.w1 * dt; s_pred.w2 += d1.w2 * dt;

    State d2 = getDerivatives(s_pred);
    s.a1 += (d1.a1 + d2.a1) * 0.5f * dt; s.a2 += (d1.a2 + d2.a2) * 0.5f * dt;
    s.w1 += (d1.w1 + d2.w1) * 0.5f * dt; s.w2 += (d1.w2 + d2.w2) * 0.5f * dt;
    return s;
}

State solveRK4(State s, float dt) {
    State k1 = getDerivatives(s);
    State s2 = s; s2.a1 += k1.a1 * dt * 0.5f; s2.a2 += k1.a2 * dt * 0.5f; s2.w1 += k1.w1 * dt * 0.5f; s2.w2 += k1.w2 * dt * 0.5f;
    State k2 = getDerivatives(s2);
    State s3 = s; s3.a1 += k2.a1 * dt * 0.5f; s3.a2 += k2.a2 * dt * 0.5f; s3.w1 += k2.w1 * dt * 0.5f; s3.w2 += k2.w2 * dt * 0.5f;
    State k3 = getDerivatives(s3);
    State s4 = s; s4.a1 += k3.a1 * dt; s4.a2 += k3.a2 * dt; s4.w1 += k3.w1 * dt; s4.w2 += k3.w2 * dt;
    State k4 = getDerivatives(s4);

    s.a1 += (k1.a1 + 2 * k2.a1 + 2 * k3.a1 + k4.a1) * dt / 6.0f;
    s.a2 += (k1.a2 + 2 * k2.a2 + 2 * k3.a2 + k4.a2) * dt / 6.0f;
    s.w1 += (k1.w1 + 2 * k2.w1 + 2 * k3.w1 + k4.w1) * dt / 6.0f;
    s.w2 += (k1.w2 + 2 * k2.w2 + 2 * k3.w2 + k4.w2) * dt / 6.0f;
    return s;
}

float calculateEnergy(const State& s) {
    float y1 = -L1 * cos(s.a1);
    float y2 = y1 - L2 * cos(s.a2);
    float V = M1 * G * y1 + M2 * G * y2;
    float T = 0.5f * M1 * (L1 * L1) * (s.w1 * s.w1) +
        0.5f * M2 * ((L1 * L1) * (s.w1 * s.w1) + (L2 * L2) * (s.w2 * s.w2) +
            2 * L1 * L2 * s.w1 * s.w2 * cos(s.a1 - s.a2));
    return T + V;
}

int main() {
    sf::ContextSettings settings;
    // В SFML 3.0 изменили название (AntiAliasingLevel -> antiAliasingLevel)
    settings.antiAliasingLevel = 8;
    auto state1 = sf::State::Windowed;

    sf::Vector2u size(1000u, 700u); 
    sf::VideoMode mode(size);

    // 3. Явно задаем стиль (заголовок + кнопка закрыть + изменить размер)
    // sf::Style::Default - это битовая маска (uint32_t)
    auto style = sf::Style::Default;

    // 4. И наконец создаем окно, передавая уже готовые переменные
    sf::RenderWindow window(sf::VideoMode({1000, 700}), "Physics Lab: Double Pendulum", style, state1, settings);
    window.setFramerateLimit(60);

    sf::Font font;
    // Если шрифта нет, программа не упадет, но текста не будет. Лучше закинь arial.ttf к .exe
    if (!font.openFromFile("C:/Windows/Fonts/arial.ttf")) {
        std::cerr << "Font not found" << std::endl;
    }

    State initialState = { 3.14159f / 2.0f, 3.14159f / 2.0f, 0, 0 };
    State state = initialState;

    float dt = 0.016f;
    float initialEnergy = calculateEnergy(state);
    int currentMethod = 3;
    std::string methodNames[] = { "", "Euler", "Trapezoidal", "Runge-Kutta 4" };

    // Храним точки хвоста в обычном векторе
    std::vector<sf::Vector2f> trailPoints;

    while (window.isOpen()) {
        // --- НОВАЯ ОБРАБОТКА СОБЫТИЙ (SFML 3.0) ---
        while (const std::optional event = window.pollEvent()) {
            if (event->is<sf::Event::Closed>()) {
                window.close();
            }
            // Проверка нажатия клавиш
            if (const auto* keyPress = event->getIf<sf::Event::KeyPressed>()) {
                if (keyPress->code == sf::Keyboard::Key::Num1) {
                    currentMethod = 1; state = initialState; trailPoints.clear(); initialEnergy = calculateEnergy(state);
                }
                if (keyPress->code == sf::Keyboard::Key::Num2) {
                    currentMethod = 2; state = initialState; trailPoints.clear(); initialEnergy = calculateEnergy(state);
                }
                if (keyPress->code == sf::Keyboard::Key::Num3) {
                    currentMethod = 3; state = initialState; trailPoints.clear(); initialEnergy = calculateEnergy(state);
                }
                if (keyPress->code == sf::Keyboard::Key::R) {
                    state = initialState; trailPoints.clear();
                }
            }
        }

        // Физика (суб-степпинг)
        int subSteps = 10;
        float subDt = dt / subSteps;
        for (int i = 0; i < subSteps; i++) {
            if (currentMethod == 1) state = solveEuler(state, subDt);
            else if (currentMethod == 2) state = solveTrapezoidal(state, subDt);
            else if (currentMethod == 3) state = solveRK4(state, subDt);
        }

        float currentEnergy = calculateEnergy(state);
        float error = std::abs(currentEnergy - initialEnergy);

        float x1 = SCALE * L1 * sin(state.a1);
        float y1 = SCALE * L1 * cos(state.a1);
        float x2 = x1 + SCALE * L2 * sin(state.a2);
        float y2 = y1 + SCALE * L2 * cos(state.a2);

        // Обновляем хвост
        trailPoints.push_back({ x2 + OFFSET_X, y2 + OFFSET_Y }); // Фигурные скобки для вектора
        if (trailPoints.size() > 1000) {
            trailPoints.erase(trailPoints.begin());
        }

        window.clear(sf::Color(30, 30, 30));

        // Отрисовка текста
        sf::Text text(font); // Новый конструктор текста
        text.setCharacterSize(18);
        text.setFillColor(sf::Color::White);

        std::stringstream ss;
        ss << "Method: " << methodNames[currentMethod] << "\n";
        ss << "Energy Error: " << std::scientific << error;
        text.setString(ss.str());
        text.setPosition({ 10.f, 10.f }); // Вектор в фигурных скобках
        window.draw(text);

        // Отрисовка хвоста
        if (trailPoints.size() > 1) {
            sf::VertexArray trace(sf::PrimitiveType::LineStrip, trailPoints.size());
            for (size_t i = 0; i < trailPoints.size(); ++i) {
                trace[i].position = trailPoints[i];
                trace[i].color = sf::Color(0, 255, 255, 150);
            }
            window.draw(trace);
        }

        // Линии маятника
        sf::Vertex line1[] = { sf::Vertex({OFFSET_X, OFFSET_Y}), sf::Vertex({x1 + OFFSET_X, y1 + OFFSET_Y}) };
        sf::Vertex line2[] = { sf::Vertex({x1 + OFFSET_X, y1 + OFFSET_Y}), sf::Vertex({x2 + OFFSET_X, y2 + OFFSET_Y}) };

        window.draw(line1, 2, sf::PrimitiveType::Lines);
        window.draw(line2, 2, sf::PrimitiveType::Lines);

        // Грузы
        sf::CircleShape m1(10); m1.setOrigin({ 10, 10 }); m1.setPosition({ x1 + OFFSET_X, y1 + OFFSET_Y }); m1.setFillColor(sf::Color::Red);
        sf::CircleShape m2(10); m2.setOrigin({ 10, 10 }); m2.setPosition({ x2 + OFFSET_X, y2 + OFFSET_Y }); m2.setFillColor(sf::Color::Red);

        window.draw(m1);
        window.draw(m2);

        window.display();
    }
    return 0;
}