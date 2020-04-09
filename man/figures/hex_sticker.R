library(hexSticker)
library(magick)

# image from : https://icon-icons.com/icon/tumor-cancer-human-biology-disease-medical-health-cell/133488
imgurl <- "man/figures/tumor.png"

sticker(imgurl, 
        package="curatedPCaData", p_size=5, p_color = "#FB62AC", p_y = 0.6,
        s_x=1, s_y=1.2, s_width=.6,
        h_fill = "white", h_color = "#2B2671",
        url = "https://github.com/Syksy/curatedPCaData", u_color = "#FB62AC", u_size = 1,
        filename="man/figures/hex.png")