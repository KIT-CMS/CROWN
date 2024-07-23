#ifndef GUARDLOGGER_H
#define GUARDLOGGER_H

#include "ROOT/RVec.hxx"
#include <Math/Vector4D.h>
#include <TObjString.h>
#include <bitset>
#include <fmt/core.h> // Include fmt library header
#include <map>
#include <mutex>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

// Specialize fmt::formatter for TString
template <> struct fmt::formatter<TString> : fmt::formatter<std::string> {
    // Use the format method of the std::string formatter
    template <typename FormatContext>
    auto format(const TString &tstring,
                FormatContext &ctx) -> decltype(ctx.out()) {
        // Directly use the formatter for std::string on TString's value
        return fmt::formatter<std::string>::format((std::string)tstring, ctx);
    }
};

// Specialize fmt::formatter for RVecs
template <typename T> struct fmt::formatter<ROOT::VecOps::RVec<T>> {
    // Parse is required by the fmt library. In this case, we don't need to
    // parse anything specific.
    constexpr auto parse(format_parse_context &ctx) -> decltype(ctx.begin()) {
        return ctx.begin();
    }

    // Implement the format method to define how to convert
    // ROOT::VecOps::RVec<T> to a string
    template <typename FormatContext>
    constexpr auto format(const ROOT::VecOps::RVec<T> &vec,
                          FormatContext &ctx) -> decltype(ctx.out()) {
        // Start with an opening bracket
        fmt::format_to(ctx.out(), "[");
        // Format each element separated by a comma
        for (size_t i = 0; i < vec.size(); ++i) {
            if (i > 0) {
                fmt::format_to(ctx.out(), ", ");
            }
            fmt::format_to(ctx.out(), "{}", vec[i]);
        }
        // End with a closing bracket
        return fmt::format_to(ctx.out(), "]\n");
    }
};

// specialization for fmt::formatter for ROOT::Math::PtEtaPhiMVector
template <> struct fmt::formatter<ROOT::Math::PtEtaPhiMVector> {
    // Parse is required by the fmt library. In this case, we don't need to
    // parse anything specific.
    constexpr auto parse(format_parse_context &ctx) -> decltype(ctx.begin()) {
        return ctx.begin();
    }

    // Implement the format method to define how to convert
    // ROOT::Math::PtEtaPhiMVector to a string
    template <typename FormatContext>
    constexpr auto format(const ROOT::Math::PtEtaPhiMVector &vec,
                          FormatContext &ctx) -> decltype(ctx.out()) {
        // Format the vector as needed, for example:
        return fmt::format_to(ctx.out(), "[Pt: {}, Eta: {}, Phi: {}, M: {}]\n",
                              vec.Pt(), vec.Eta(), vec.Phi(), vec.M());
    }
};

// Generic formatter function for std::bitset of any size
template <std::size_t N> struct BitsetFormatter {
    static auto format(fmt::format_context &ctx,
                       const std::bitset<N> &b) -> decltype(ctx.out()) {
        fmt::format_to(ctx.out(), "[");
        for (size_t i = 0; i < b.size(); ++i) {
            fmt::format_to(ctx.out(), "{}", b.test(i) ? '1' : '0');
        }
        return fmt::format_to(ctx.out(), "]\n");
    }
};

// Register the formatter with fmt library for any size of std::bitset
template <std::size_t N>
struct fmt::formatter<std::bitset<N>> : formatter<std::string> {
    template <typename FormatContext>
    auto format(const std::bitset<N> &b,
                FormatContext &ctx) -> decltype(ctx.out()) {
        return BitsetFormatter<N>::format(ctx, b);
    }
};

class Logger {
  public:
    static std::shared_ptr<spdlog::logger> get(std::string name) {
        if (getInstance()._loggers.count(name) == 0) {
            std::vector<spdlog::sink_ptr> sinkVector;
            sinkVector.push_back(
                std::make_shared<spdlog::sinks::stdout_color_sink_st>());

            // check if file logging is enabled
            if (getInstance()._fileName)
                sinkVector.push_back(
                    std::make_shared<spdlog::sinks::basic_file_sink_st>(
                        *getInstance()._fileName));

            auto newLogger = std::make_shared<spdlog::logger>(
                name, begin(sinkVector), end(sinkVector));
            newLogger->set_level(convertLevelToSpdlog(getInstance()._level));
            getInstance()._loggers[name] = newLogger;
        }

        return getInstance()._loggers[name];
    }
    enum class LogLevel { DEBUG, INFO, WARN, ERR, CRITICAL, OFF };
    static void setLevel(LogLevel level) {
        getInstance()._level = level;

        // set level globally (probably superfluous..)
        spdlog::set_level(convertLevelToSpdlog(level));

        // set level for all active loggers
        for (auto &[key, logger] : getInstance()._loggers)
            logger->set_level(convertLevelToSpdlog(level));
    }
    static void enableFileLogging(std::string filename) {
        getInstance()._fileName = std::make_unique<std::string>(filename);
        for (auto &[key, logger] : getInstance()._loggers) {
            // if there is less than two sinks, add a file sink
            if (logger->sinks().size() < 2)
                logger->sinks().push_back(
                    std::make_shared<spdlog::sinks::basic_file_sink_st>(
                        *getInstance()._fileName));
        }
    }
    ~Logger() {
        _loggers.clear();  // Clear the map of loggers
        _fileName.reset(); // Reset the unique_ptr, effectively deleting the
                           // managed object
    }

  private:
    static Logger &getInstance() {
        static Logger instance;
        return instance;
    }
    LogLevel _level{LogLevel::INFO};
    static spdlog::level::level_enum convertLevelToSpdlog(LogLevel level) {
        switch (level) {
        case LogLevel::DEBUG:
            return spdlog::level::debug;
        case LogLevel::INFO:
            return spdlog::level::info;
        case LogLevel::WARN:
            return spdlog::level::warn;
        case LogLevel::ERR:
            return spdlog::level::err;
        case LogLevel::CRITICAL:
            return spdlog::level::critical;
        case LogLevel::OFF:
            return spdlog::level::off;
        default:
            return spdlog::level::info;
        }
    }
    std::unique_ptr<std::string> _fileName{};
    std::map<std::string, std::shared_ptr<spdlog::logger>> _loggers;
};

#endif /* GUARDLOGGER_H */
