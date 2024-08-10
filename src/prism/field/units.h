#pragma once

#include <cstdint>

namespace prism::field::units {
struct SystemOfUnits {
  public:
    void inline setTime(int8_t power) noexcept { _time = power; }
    void inline setLength(int8_t power) noexcept { _length = power; }
    void inline setMass(int8_t power) noexcept { _mass = power; }
    void inline setCurrent(int8_t power) noexcept { _current = power; }
    void inline setTemperature(int8_t power) noexcept { _temperature = power; }
    void inline setSubstanceAmount(int8_t power) noexcept { _substance_amount = power; }
    void inline setLuminousIntensity(int8_t power) noexcept { _luminous_intensity = power; }

    auto time() const noexcept -> int8_t { return _time; }
    auto length() const noexcept -> int8_t { return _length; }
    auto mass() const noexcept -> int8_t { return _mass; }
    auto current() const noexcept -> int8_t { return _current; }
    auto temperature() const noexcept -> int8_t { return _temperature; }
    auto substanceAmount() const noexcept -> int8_t { return _substance_amount; }
    auto luminousIntensity() const noexcept -> int8_t { return _luminous_intensity; }

  private:
    int8_t _time {0};
    int8_t _length {0};
    int8_t _mass {0};
    int8_t _current {0};
    int8_t _temperature {0};
    int8_t _substance_amount {0};
    int8_t _luminous_intensity {0};
};

class Measurable {
  public:
    Measurable() = default;
    Measurable(const Measurable&) noexcept = default;
    Measurable(Measurable&&) noexcept = default;
    auto operator=(const Measurable&) -> Measurable& = default;
    auto operator=(Measurable&&) noexcept -> Measurable& = default;
    virtual ~Measurable() = default;

    virtual auto units() -> SystemOfUnits& { return _units; }
    virtual auto units() const -> const SystemOfUnits& { return _units; }

  private:
    SystemOfUnits _units;
};

class VelocityUnit : public Measurable {
  public:
    VelocityUnit() {
        units().setLength(1);
        units().setTime(-1);
    }
};

class PressureUnit : public Measurable {
  public:
    PressureUnit() {
        units().setMass(1);
        units().setLength(1);
        units().setTime(-2);
    }
};

using Pascal = PressureUnit;

} // namespace prism::field::units