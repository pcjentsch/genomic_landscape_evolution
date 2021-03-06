

#estimated from Wilks, S. H., Mühlemann, B., Shen, X., Türeli, S., LeGresley, E. B., Netzl, A., Caniza, M. A., Chacaltana-Huarcaya, J. N., Daniell, X., Datto, M. B., Denny, T. N., Drosten, C., Fouchier, R. A. M., Garcia, P. J., Halfmann, P. J., Jassem, A., Jones, T. C., Kawaoka, Y., Krammer, F., … Smith, D. J. (2022). Mapping SARS-CoV-2 antigenic relationships and serological responses [Preprint]. Immunology. https://doi.org/10.1101/2022.01.28.477987
#figure 2
const antigenic_map = [
    ("B.1.1.28.1", (4.2, 5.1, 0.4),),
    ("B.1.1.7", (2.7, 3.8, 0.4),),
    ("B.1", (2.7, 3.2, 0.4),), #D614G, calling this B.1 for now 
    ("B.1.429", (3.4, 2.2, 0.4),),
    ("B.1.617.2", (2.25, 1.75, 0.4),),
    ("B.1.1.1.37", (4.75, 2.25, 0.4),),
    ("B.1.617.1", (5.6, 2.1, 0.4),),
    ("B.1.1.529", (7.8, 1.5, 0.4),),
    ("B.1.621", (5.8, 4.7, 0.4),),
    ("B.1.351", (5.3, 5.6, 0.4),),
]

const antigenic_map_paper = Dict(
    missing => missing,
    "B.1.1.28.1" => (4.2, 5.1, 0.4),
    "B.1.1.7" => (2.7, 3.8, 0.4),
    "D614G" => (2.7, 3.2, 0.4),
    "B.1.429" => (3.4, 2.2, 0.4),
    "B.1.617.2" => (2.25, 1.75, 0.4),
    "B.1.1.1.37" => (4.75, 2.25, 0.4),
    "B.1.617.1" => (5.6, 2.1, 0.4),
    "B.1.1.529" => (7.8, 1.5, 0.4),
    "B.1.621" => (5.8, 4.7, 0.4),
    "B.1.351" => (5.3, 5.6, 0.4),
    "B.1.526+E484K" => (5.6, 3.4, 0.4),
    "B.1.617.2(AY.2)+K417N" => (2.5, 1.75, 0.2),
    "B.1.617.2(AY.1)+K417N" => (2.2, 1.75, 0.2),
    "B.1.617.2+K417N" => (2.3, 2.0, 0.2),
    "B.1.617.2(AY.3)+E484Q" => (5.5, 1.8, 0.2),
    "B.1.1.7+E484K" => (5.2, 4.1, 0.2),
)

