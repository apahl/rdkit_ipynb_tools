function toggle_clusters() {
    if (window.flag_clusters_expanded) {
        $('ul ul').hide();
        window.flag_clusters_expanded = false;
    } else {
        $('ul ul').show();
        window.flag_clusters_expanded = true;
    }
}

function toggle_histograms() {
    if (window.flag_histograms_shown) {
        $('.histogram').hide();
        window.flag_histograms_shown = false;
    } else {
        $('.histogram').show();
        window.flag_histograms_shown = true;
    }
}
