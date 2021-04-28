namespace metfilter {

auto ApplyMetFilter(auto df, const std::string &flagname,
                    const std::string &filtername) {
    return df.Filter([](const bool flag) { return flag; }, {flagname},
                     filtername);
}

} // namespace metfilter