#include <SFML/Graphics.hpp>
#include <cmath>
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>
#include <optional>
#include <iostream>
#include <cstdint>

// I tryed to discribe all parts of code in comments, if you have any questions about it, please ask me. I will be happy to explain it to you)
// p.s sorry for russian comments, it was more understandable for my team and our's presentation. If you want, I can translate them to English as well)

// SIMULATED "BUILT-IN" LIBRARY SECTION
// Since C++ STL doesn't have a built-in ODE solver, we encapsulate
// the physics logic in a namespace to simulate a library structure.
namespace PhysicsEngine {

    // Physical constants and simulation parameters
    struct Params {
        float G = 9.81f;       // Gravity
        float M1 = 1.0f;       // Mass of first pendulum
        float M2 = 1.0f;       // Mass of second pendulum
        float L1 = 1.0f;       // Length of first rod
        float L2 = 1.0f;       // Length of second rod
        float timeScale = 1.0f;// Simulation speed factor
    };

    // System state: angles and angular velocities
    struct State {
        float a1 = 0; float a2 = 0; // Angles
        float w1 = 0; float w2 = 0; // Angular Velocities
    };

    // INTERNAL HELPER FUNCTION
    // Calculates derivatives based on Lagrangian mechanics.
    State getDerivatives(const State& s, const Params& p) {
        State d;
        // Связь между углом и скоростью.
        // Производная от угла (позиции) — это всегда угловая скорость.
        d.a1 = s.w1;
        d.a2 = s.w2;

        // Строка 7: Разница углов.
        // От нее зависит, как сильно маятники тянут друг друга.
        // Если del = 0, они висят в линию. Если del = 90 градусов, сила передачи меняется.
        float del = s.a1 - s.a2;

        // Общий знаменатель.
        // Эта формула получена из решения системы линейных уравнений Лагранжа.
        float commonDenom = 2 * p.M1 + p.M2 - p.M2 * cos(2 * del);

        // Числитель для ПЕРВОГО маятника(num1).
        // Это сумма всех сил, действующих на верхний стержень:
        // 1. -p.G * ... * sin(s.a1) -> Гравитация тянет его вниз.
        // 2. Второе слагаемое -> Влияние гравитации второго шара (через рычаг).
        // 3. Третье слагаемое (с s.w1, s.w2) -> Центробежные силы и инерция.
        float num1 = -p.G * (2 * p.M1 + p.M2) * sin(s.a1)
            - p.M2 * p.G * sin(s.a1 - 2 * s.a2)
            - 2 * sin(del) * p.M2 * (s.w2 * s.w2 * p.L2 + s.w1 * s.w1 * p.L1 * cos(del));

        // Числитель для ВТОРОГО маятника (num2).
        // Силы, действующие на нижний стержень.
        float num2 = 2 * sin(del) * (s.w1 * s.w1 * p.L1 * (p.M1 + p.M2) + p.G * (p.M1 + p.M2) * cos(s.a1)
            + s.w2 * s.w2 * p.L2 * p.M2 * cos(del));

        // Финальный расчет УГЛОВОГО УСКОРЕНИЯ.
        // Ускорение = Сила / Массу (грубо говоря, по 2 закону Ньютона).
        // Здесь мы делим числители на знаменатель и длину стержня.
        d.w1 = num1 / (p.L1 * commonDenom);
        d.w2 = num2 / (p.L2 * commonDenom);
        return d;
    }

    // TRAPEZOIDAL RULE (Heun's Method)
    // Less accurate, good for demonstration.
    State solveTrapezoidal(State s, float dt, const Params& p) {
        // Шаг 1: ПРЕДИКТОР (Предсказание по Эйлеру).
        // Мы спрашиваем: "Какой наклон (k1) прямо сейчас?"
        State k1 = getDerivatives(s, p);

        State s_guess = s;
        s_guess.a1 += k1.a1 * dt; s_guess.a2 += k1.a2 * dt;
        s_guess.w1 += k1.w1 * dt; s_guess.w2 += k1.w2 * dt;

        // Мы делаем "черновой" шаг вперед, используя текущий наклон.
        // s_guess — это то, где мы окажемся, если пойдем тупо по прямой.
        State k2 = getDerivatives(s_guess, p);

		// Вычисляем средний наклон между началом и концом отрезка
        s.a1 += (k1.a1 + k2.a1) * 0.5f * dt;
        s.a2 += (k1.a2 + k2.a2) * 0.5f * dt;
        s.w1 += (k1.w1 + k2.w1) * 0.5f * dt;
        s.w2 += (k1.w2 + k2.w2) * 0.5f * dt;

        return s;
    }

    // SOLVER 2: RUNGE-KUTTA 4 (RK4)
    // This is the industry standard for physics simulations.
    // Replaces Simpson's rule (which is for integrals) in the context of ODEs.
    State solveRK4(State s, float dt, const Params& p) {
        // Делаем первый расчет с инзачального места
        State k1 = getDerivatives(s, p);

        // делаем пол-шага вперед относительного прошлого и смотрим
        State s2 = s;
        s2.a1 += k1.a1 * dt * 0.5f; s2.a2 += k1.a2 * dt * 0.5f;
        s2.w1 += k1.w1 * dt * 0.5f; s2.w2 += k1.w2 * dt * 0.5f;
        State k2 = getDerivatives(s2, p);

        // делаем еще пол шага вперед, опять же относительно прошлого
        State s3 = s;
        s3.a1 += k2.a1 * dt * 0.5f; s3.a2 += k2.a2 * dt * 0.5f;
        s3.w1 += k2.w1 * dt * 0.5f; s3.w2 += k2.w2 * dt * 0.5f;
        State k3 = getDerivatives(s3, p);

        // делаем один шаг a to b (из прошлого кадра в новый), выглядит как отрезок
        State s4 = s;
        s4.a1 += k3.a1 * dt; s4.a2 += k3.a2 * dt;
        s4.w1 += k3.w1 * dt; s4.w2 += k3.w2 * dt;
        State k4 = getDerivatives(s4, p);

        //  тут же мы берем среднее из них, но увеличиваем вес второго и третьего в 2 раза т.к они показывают угол близкий к реальности 
        s.a1 += (k1.a1 + 2 * k2.a1 + 2 * k3.a1 + k4.a1) * dt / 6.0f;
        s.a2 += (k1.a2 + 2 * k2.a2 + 2 * k3.a2 + k4.a2) * dt / 6.0f;
        s.w1 += (k1.w1 + 2 * k2.w1 + 2 * k3.w1 + k4.w1) * dt / 6.0f;
        s.w2 += (k1.w2 + 2 * k2.w2 + 2 * k3.w2 + k4.w2) * dt / 6.0f;

        return s;
    }

    // Calculates Total Mechanical Energy (T + V) for Error Analysis
    float getEnergy(const State& s, const Params& p) {
        // Potential Energy (V = mgh)
        float y1 = -p.L1 * cos(s.a1);
        float y2 = y1 - p.L2 * cos(s.a2);
		// mgh для первого маятника + mgh для второго маятника
        float V = p.M1 * p.G * y1 + p.M2 * p.G * y2;

		// Kinetic Energy (T = 0.5 * m * v^2)
        float term1 = 0.5f * p.M1 * (p.L1 * p.L1) * (s.w1 * s.w1);

        // Для второго маятника все СЛОЖНО.
        // Его скорость зависит от движения первого + его собственного вращения.
        // Тут используется теорема косинусов для сложения векторов скоростей.
        // Формула: 0.5 * m2 * (v1^2 + v2^2 + 2*v1*v2*cos(a1-a2))
        float term2 = 0.5f * p.M2 * ((p.L1 * p.L1) * (s.w1 * s.w1) +
            (p.L2 * p.L2) * (s.w2 * s.w2) +
            2 * p.L1 * p.L2 * s.w1 * s.w2 * cos(s.a1 - s.a2));

        float T = term1 + term2;
        return T + V;
    }
}


int main() {
    // SFML 3.0 Window Setup
    sf::ContextSettings settings;
    settings.antiAliasingLevel = 8;
    sf::Vector2u windowSize(1200u, 800u);
    sf::VideoMode mode(windowSize);
    sf::String title("Double Pendulum Lab: RK4 vs Trapezoidal");
    auto style = static_cast<uint32_t>(sf::Style::Default);
    auto windowState = sf::State::Windowed;
    sf::RenderWindow window(mode, title, style, windowState, settings);
    window.setFramerateLimit(60);

    // Load font
    sf::Font font;
    if (!font.openFromFile("C:/Windows/Fonts/arial.ttf")) {
        std::cerr << "Arial font not found! Text will not be displayed." << std::endl;
    }

    PhysicsEngine::Params params;
    PhysicsEngine::State initialState = { 3.14159f / 2.0f, 3.14159f / 2.0f, 0, 0 };
    PhysicsEngine::State state = initialState;

    // Track initial energy
    float initialEnergy = PhysicsEngine::getEnergy(state, params);

    std::vector<sf::Vector2f> trailPoints;
    const float SCALE = 150.0f;
    const sf::Vector2f OFFSET(windowSize.x / 2.0f, windowSize.y / 3.0f);

    int selectedItem = 0;
    const int ITEM_COUNT = 6;
    std::string paramNames[] = { "Gravity (G)", "Mass 1 (kg)", "Mass 2 (kg)", "Length 1 (m)", "Length 2 (m)", "Time Scale" };
    bool isPaused = true;

    // Toggle for solving method: true = RK4, false = Trapezoidal
    bool useRK4 = true;

    while (window.isOpen()) {
        // Event Polling
        while (const std::optional event = window.pollEvent()) {
            if (event->is<sf::Event::Closed>()) window.close();

            if (const auto* keyPress = event->getIf<sf::Event::KeyPressed>()) {
                if (keyPress->code == sf::Keyboard::Key::R) {
                    state = initialState;
                    trailPoints.clear();
                    initialEnergy = PhysicsEngine::getEnergy(state, params);
                }
                if (keyPress->code == sf::Keyboard::Key::P) isPaused = !isPaused;

                // Toggle Method with 'M' key
                if (keyPress->code == sf::Keyboard::Key::M) {
                    useRK4 = !useRK4;
                    // Reset to clearly see the difference in stability
                    state = initialState;
                    trailPoints.clear();
                    initialEnergy = PhysicsEngine::getEnergy(state, params);
                }

                if (keyPress->code == sf::Keyboard::Key::Up) {
                    selectedItem--; if (selectedItem < 0) selectedItem = ITEM_COUNT - 1;
                }
                if (keyPress->code == sf::Keyboard::Key::Down) {
                    selectedItem++; if (selectedItem >= ITEM_COUNT) selectedItem = 0;
                }
            }
        }

        // Parameter control
        float changeSpeed = 0.0f;
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Left)) changeSpeed = -0.05f;
        else if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Right)) changeSpeed = 0.05f;

        if (changeSpeed != 0.0f) {
            if (selectedItem == 0) params.G += changeSpeed;
            if (selectedItem == 1) { params.M1 += changeSpeed; if (params.M1 < 0.1f) params.M1 = 0.1f; }
            if (selectedItem == 2) { params.M2 += changeSpeed; if (params.M2 < 0.1f) params.M2 = 0.1f; }
            if (selectedItem == 3) { params.L1 += changeSpeed * 0.1f; if (params.L1 < 0.1f) params.L1 = 0.1f; }
            if (selectedItem == 4) { params.L2 += changeSpeed * 0.1f; if (params.L2 < 0.1f) params.L2 = 0.1f; }
            if (selectedItem == 5) { params.timeScale += changeSpeed * 0.1f; if (params.timeScale < 0.0f) params.timeScale = 0.0f; }

            if (selectedItem == 3 || selectedItem == 4) trailPoints.clear();
            initialEnergy = PhysicsEngine::getEnergy(state, params);
        }

        if (!isPaused) {
            float dt = 0.016f * params.timeScale;
            int subSteps = 10;
            float subDt = dt / subSteps;

            for (int i = 0; i < subSteps; i++) {
                if (useRK4) {
                    state = PhysicsEngine::solveRK4(state, subDt, params);
                }
                else {
                    state = PhysicsEngine::solveTrapezoidal(state, subDt, params);
                }
            }

            // Trail update
            float x1 = SCALE * params.L1 * sin(state.a1);
            float x2 = x1 + SCALE * params.L2 * sin(state.a2);
            float y1 = SCALE * params.L1 * cos(state.a1);
            float y2 = y1 + SCALE * params.L2 * cos(state.a2);

            trailPoints.push_back({ x2 + OFFSET.x, y2 + OFFSET.y });
            if (trailPoints.size() > 1000) trailPoints.erase(trailPoints.begin());
        }

        // --- RENDER ---
        window.clear(sf::Color(20, 20, 25));

        // Draw Trail
        if (trailPoints.size() > 1) {
            sf::VertexArray trace(sf::PrimitiveType::LineStrip, trailPoints.size());
            for (size_t i = 0; i < trailPoints.size(); ++i) {
                trace[i].position = trailPoints[i];
                float alpha = 255.0f * ((float)i / trailPoints.size());
                trace[i].color = sf::Color(0, 255, 255, static_cast<uint8_t>(alpha));
            }
            window.draw(trace);
        }

        // Calculate visual positions
        float x1 = SCALE * params.L1 * sin(state.a1);
        float y1 = SCALE * params.L1 * cos(state.a1);
        float x2 = x1 + SCALE * params.L2 * sin(state.a2);
        float y2 = y1 + SCALE * params.L2 * cos(state.a2);

        sf::Vector2f center = OFFSET;
        sf::Vector2f pos1 = { x1 + OFFSET.x, y1 + OFFSET.y };
        sf::Vector2f pos2 = { x2 + OFFSET.x, y2 + OFFSET.y };

        // Draw Rods
        sf::Vertex line1[] = { sf::Vertex(center, sf::Color::White), sf::Vertex(pos1, sf::Color::White) };
        sf::Vertex line2[] = { sf::Vertex(pos1, sf::Color::White), sf::Vertex(pos2, sf::Color::White) };
        window.draw(line1, 2, sf::PrimitiveType::Lines);
        window.draw(line2, 2, sf::PrimitiveType::Lines);

        // Draw Masses
        float r1 = 10.0f + params.M1 * 2.0f;
        float r2 = 10.0f + params.M2 * 2.0f;
        sf::CircleShape m1(r1); m1.setOrigin({ r1, r1 }); m1.setPosition(pos1); m1.setFillColor(sf::Color::Red);
        sf::CircleShape m2(r2); m2.setOrigin({ r2, r2 }); m2.setPosition(pos2); m2.setFillColor(sf::Color::Red);
        window.draw(m1);
        window.draw(m2);

        sf::CircleShape hub(5); hub.setOrigin({ 5, 5 }); hub.setPosition(center); hub.setFillColor(sf::Color::White);
        window.draw(hub);

        // --- GUI ---
        sf::RectangleShape panel({ 320.0f, 800.0f });
        panel.setFillColor(sf::Color(0, 0, 0, 180));
        window.draw(panel);

        sf::Text guiText(font);
        guiText.setCharacterSize(18);
        guiText.setString("CONTROLS:\n[Up/Down] Select\n[Left/Right] Change\n[M] Switch Method\n[P] Pause  [R] Reset\n------------------");
        guiText.setPosition({ 10.f, 10.f });
        guiText.setFillColor(sf::Color::White);
        window.draw(guiText);

        float startY = 150.0f;
        float stepY = 35.0f;
        float values[] = { params.G, params.M1, params.M2, params.L1, params.L2, params.timeScale };

        for (int i = 0; i < ITEM_COUNT; i++) {
            std::stringstream ss;
            if (i == selectedItem) ss << "> "; else ss << "  ";
            ss << paramNames[i] << ": " << std::fixed << std::setprecision(2) << values[i];
            guiText.setString(ss.str());
            guiText.setPosition({ 10.f, startY + i * stepY });
            if (i == selectedItem) guiText.setFillColor(sf::Color::Yellow);
            else guiText.setFillColor(sf::Color(180, 180, 180));
            window.draw(guiText);
        }

        // --- ANALYSIS ---
        float currentEnergy = PhysicsEngine::getEnergy(state, params);
        float error = std::abs(currentEnergy - initialEnergy);

        float analysisY = startY + ITEM_COUNT * stepY + 30.0f;

        // Show which method is active
        std::string methodStr = useRK4 ? "RK4" : "Trapezoidal";
        sf::Text analysisTitle(font, " METHOD: " + methodStr + " ---", 18);
        analysisTitle.setFillColor(useRK4 ? sf::Color::Cyan : sf::Color::Magenta);
        analysisTitle.setPosition({ 10.f, analysisY });
        window.draw(analysisTitle);

        std::stringstream ssAnalysis;
        ssAnalysis << "Energy: " << std::fixed << std::setprecision(3) << currentEnergy << " J\n";
        ssAnalysis << "Error:  " << std::scientific << std::setprecision(2) << error << "\n";

        // Status logic
        std::string status;
        sf::Color statusColor;

        if (error < 0.1f) {
            status = "STABLE";
            statusColor = sf::Color::Green;
        }
        else {
            status = "DRIFT / UNSTABLE";
            statusColor = sf::Color::Red;
        }

        ssAnalysis << "Status: " << status;

        sf::Text analysisText(font, ssAnalysis.str(), 16);
        analysisText.setFillColor(statusColor);
        analysisText.setPosition({ 10.f, analysisY + 30.0f });
        window.draw(analysisText);

        if (isPaused) {
            sf::Text pauseText(font, "PAUSED", 40);
            pauseText.setFillColor(sf::Color::Red);
            pauseText.setPosition({ 600.f, 50.f });
            window.draw(pauseText);
        }

        window.display();
    }
    return 0;
}