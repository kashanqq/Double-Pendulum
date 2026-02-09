#include <SFML/Graphics.hpp>
#include <cmath>
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>
#include <optional>
#include <iostream>

struct PhysicsParams {
    float G = 9.81f;
    float M1 = 1.0f;       
    float M2 = 1.0f;      
    float L1 = 1.0f;       
    float L2 = 1.0f;       
    float timeScale = 1.0f;
};

// State of the system: angles (a1, a2) and angular velocities (w1, w2)
struct State {
    float a1 = 0; float a2 = 0;
    float w1 = 0; float w2 = 0;
};

// Formusla for derivatives based on the equations of motion for a double pendulum
State getDerivatives(const State& s, const PhysicsParams& p) {
    State d;
    d.a1 = s.w1;
    d.a2 = s.w2;

    float del = s.a1 - s.a2;
    float den1 = (p.M1 + p.M2) * p.L1 - p.M2 * p.L1 * cos(del) * cos(del);
    float den2 = (p.L2 / p.L1) * den1;

    float num1 = -p.G * (2 * p.M1 + p.M2) * sin(s.a1)
        - p.M2 * p.G * sin(s.a1 - 2 * s.a2)
        - 2 * sin(del) * p.M2 * (s.w2 * s.w2 * p.L2 + s.w1 * s.w1 * p.L1 * cos(del));

    float num2 = 2 * sin(del) * (s.w1 * s.w1 * p.L1 * (p.M1 + p.M2) + p.G * (p.M1 + p.M2) * cos(s.a1)
        + s.w2 * s.w2 * p.L2 * p.M2 * cos(del));

    d.w1 = num1 / den1;
    d.w2 = num2 / den2;
    return d;
}

// RK4 Integration method to update the state of the system
State solveRK4(State s, float dt, const PhysicsParams& p) {
    State k1 = getDerivatives(s, p);
    State s2 = s; s2.a1 += k1.a1 * dt * 0.5f; s2.a2 += k1.a2 * dt * 0.5f; s2.w1 += k1.w1 * dt * 0.5f; s2.w2 += k1.w2 * dt * 0.5f;
    State k2 = getDerivatives(s2, p);
    State s3 = s; s3.a1 += k2.a1 * dt * 0.5f; s3.a2 += k2.a2 * dt * 0.5f; s3.w1 += k2.w1 * dt * 0.5f; s3.w2 += k2.w2 * dt * 0.5f;
    State k3 = getDerivatives(s3, p);
    State s4 = s; s4.a1 += k3.a1 * dt; s4.a2 += k3.a2 * dt; s4.w1 += k3.w1 * dt; s4.w2 += k3.w2 * dt;
    State k4 = getDerivatives(s4, p);

    s.a1 += (k1.a1 + 2 * k2.a1 + 2 * k3.a1 + k4.a1) * dt / 6.0f;
    s.a2 += (k1.a2 + 2 * k2.a2 + 2 * k3.a2 + k4.a2) * dt / 6.0f;
    s.w1 += (k1.w1 + 2 * k2.w1 + 2 * k3.w1 + k4.w1) * dt / 6.0f;
    s.w2 += (k1.w2 + 2 * k2.w2 + 2 * k3.w2 + k4.w2) * dt / 6.0f;
    return s;
}

int main() {
    
    sf::ContextSettings settings;
    settings.antiAliasingLevel = 8;
	sf::Vector2u windowSize(1200u, 800u);
    sf::VideoMode mode(windowSize);
    sf::String title("Chaos Engine: Interactive Lab");
    auto style = static_cast<uint32_t>(sf::Style::Default);
    auto windowState = sf::State::Windowed;
    sf::RenderWindow window(mode, title, style, windowState, settings);
    window.setFramerateLimit(60);

    sf::Font font;
    if (!font.openFromFile("C:/Windows/Fonts/arial.ttf")) {
        std::cerr << "ERROR: Arial font not found!" << std::endl;
    }

    // Variables
    PhysicsParams params; 

    // Start position (changeable)
    State initialState = { 3.14159f, 3.14159f, 0, 0 };
    State state = initialState;

    std::vector<sf::Vector2f> trailPoints;
	const float SCALE = 150.0f; // Scale factor for rendering (pixels per meter)
	const sf::Vector2f OFFSET(windowSize.x / 2.0f, windowSize.y / 2.0f); // The center of the pendulum in window coordinates
	std::cout << OFFSET.x << " " << OFFSET.y << std::endl;

	// Interface variables
    int selectedItem = 0;
    const int ITEM_COUNT = 6;
	// Names for the parameters in the GUI
    std::string paramNames[] = { "Gravity (G)", "Mass 1 (kg)", "Mass 2 (kg)", "Length 1 (m)", "Length 2 (m)", "Time Scale" };

    bool isPaused = false;

    while (window.isOpen()) {
		// Start Input Handling
        while (const std::optional event = window.pollEvent()) {
            if (event->is<sf::Event::Closed>()) window.close();

            if (const auto* keyPress = event->getIf<sf::Event::KeyPressed>()) {
                // R - Restart
                if (keyPress->code == sf::Keyboard::Key::R) {
                    state = initialState;
                    trailPoints.clear();
                }
				// P - Pause/Unpause
                if (keyPress->code == sf::Keyboard::Key::P) {
                    isPaused = !isPaused;
                }
				// Having Up/Down keys to navigate the menu
                if (keyPress->code == sf::Keyboard::Key::Up) {
                    selectedItem--;
                    if (selectedItem < 0) selectedItem = ITEM_COUNT - 1;
                }
                if (keyPress->code == sf::Keyboard::Key::Down) {
                    selectedItem++;
                    if (selectedItem >= ITEM_COUNT) selectedItem = 0;
                }
            }
        }
		// End Input Handling

		float changeSpeed = 0.05f; // Step for changing parameters
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Left)) changeSpeed = -0.05f;
        else if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Right)) changeSpeed = 0.05f;
        else changeSpeed = 0.0f;

        if (changeSpeed != 0.0f) {
            if (selectedItem == 0) params.G += changeSpeed;
            if (selectedItem == 1) { params.M1 += changeSpeed; if (params.M1 < 0.1f) params.M1 = 0.1f; }
            if (selectedItem == 2) { params.M2 += changeSpeed; if (params.M2 < 0.1f) params.M2 = 0.1f; }
            if (selectedItem == 3) { params.L1 += changeSpeed * 0.1f; if (params.L1 < 0.1f) params.L1 = 0.1f; }
            if (selectedItem == 4) { params.L2 += changeSpeed * 0.1f; if (params.L2 < 0.1f) params.L2 = 0.1f; }
            if (selectedItem == 5) { params.timeScale += changeSpeed * 0.1f; if (params.timeScale < 0.0f) params.timeScale = 0.0f; }

            if (selectedItem == 3 || selectedItem == 4) trailPoints.clear();
        }

		// Update physics only if not paused
        if (!isPaused) {
            float dt = 0.016f * params.timeScale;
            int subSteps = 10;
            float subDt = dt / subSteps;

            for (int i = 0; i < subSteps; i++) {
                state = solveRK4(state, subDt, params);
            }

			// Logic for trail: we calculate the position of the second mass and add it to the trail points. We also limit the number of points to avoid memory issues
            float x1 = SCALE * params.L1 * sin(state.a1);
            float y1 = SCALE * params.L1 * cos(state.a1);
            float x2 = x1 + SCALE * params.L2 * sin(state.a2);
            float y2 = y1 + SCALE * params.L2 * cos(state.a2);

            trailPoints.push_back({ x2 + OFFSET.x, y2 + OFFSET.y });
            if (trailPoints.size() > 1500) trailPoints.erase(trailPoints.begin());
        }

		// Start Drawing
        window.clear(sf::Color(20, 20, 25)); 

		// Draw the trail of the second mass. We use a gradient color to make it visually appealing, fading out older points
        if (trailPoints.size() > 1) {
            sf::VertexArray trace(sf::PrimitiveType::LineStrip, trailPoints.size());
            for (size_t i = 0; i < trailPoints.size(); ++i) {
                trace[i].position = trailPoints[i];
                // √радиент цвета от прозрачного к €ркому
                float alpha = 255.0f * ((float)i / trailPoints.size());
                trace[i].color = sf::Color(0, 255, 255, static_cast<uint8_t>(alpha));
            }
            window.draw(trace);
        }

		// Calculate positions of masses
        float x1 = SCALE * params.L1 * sin(state.a1);
        float y1 = SCALE * params.L1 * cos(state.a1);
        float x2 = x1 + SCALE * params.L2 * sin(state.a2);
        float y2 = y1 + SCALE * params.L2 * cos(state.a2);

        sf::Vector2f center = OFFSET;
        sf::Vector2f pos1 = { x1 + OFFSET.x, y1 + OFFSET.y };
        sf::Vector2f pos2 = { x2 + OFFSET.x, y2 + OFFSET.y };

		// Kernels (Lines)
        sf::Vertex line1[] = { sf::Vertex(center, sf::Color::White), sf::Vertex(pos1, sf::Color::White) };
        sf::Vertex line2[] = { sf::Vertex(pos1, sf::Color::White), sf::Vertex(pos2, sf::Color::White) };
        window.draw(line1, 2, sf::PrimitiveType::Lines);
        window.draw(line2, 2, sf::PrimitiveType::Lines);

		// Masses (Circles) - depending on mass value, we can change the radius for better visualization
        float r1 = 10.0f + params.M1 * 2.0f;
        float r2 = 10.0f + params.M2 * 2.0f;

        sf::CircleShape m1(r1); m1.setOrigin({ r1, r1 }); m1.setPosition(pos1); m1.setFillColor(sf::Color::Red);
        sf::CircleShape m2(r2); m2.setOrigin({ r2, r2 }); m2.setPosition(pos2); m2.setFillColor(sf::Color::Red);
        window.draw(m1);
        window.draw(m2);

		// The hub (center point)
        sf::CircleShape hub(5); hub.setOrigin({ 5, 5 }); hub.setPosition(center); hub.setFillColor(sf::Color::White);
        window.draw(hub);


		// Start GUI
        sf::RectangleShape panel({ 300.0f, 800.0f });
        panel.setFillColor(sf::Color(0, 0, 0, 150));
        window.draw(panel);

        sf::Text guiText(font);
        guiText.setCharacterSize(20);

		// Description text
        guiText.setString("CONTROLS:\n[Up/Down] Select\n[Left/Right] Change\n[P] Pause\n[R] Reset\n------------------");
        guiText.setPosition({ 10.f, 10.f });
        guiText.setFillColor(sf::Color::White);
        window.draw(guiText);

        float startY = 150.0f;
        float stepY = 40.0f;

        float values[] = { params.G, params.M1, params.M2, params.L1, params.L2, params.timeScale };

        for (int i = 0; i < ITEM_COUNT; i++) {
            std::stringstream ss;

            if (i == selectedItem) ss << "> ";
            else ss << "  ";

            ss << paramNames[i] << ": " << std::fixed << std::setprecision(2) << values[i];

            guiText.setString(ss.str());
            guiText.setPosition({ 10.f, startY + i * stepY });

            if (i == selectedItem) guiText.setFillColor(sf::Color::Yellow);
            else guiText.setFillColor(sf::Color(180, 180, 180));

            window.draw(guiText);
        }

        if (isPaused) {
            sf::Text pauseText(font, "PAUSED", 40);
            pauseText.setFillColor(sf::Color::Red);
            pauseText.setPosition({ 550.f, 50.f });
            window.draw(pauseText);
        }

        window.display();
        // End GUI
    }
    return 0;
}