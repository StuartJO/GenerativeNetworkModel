function rgb_a = make_alpha_rgb(RGB,alpha)

rgb_adjust = 1-RGB;

rgb_a = (rgb_adjust.*alpha)+RGB;