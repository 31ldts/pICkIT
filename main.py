from src.pickit.analyze_interactions import AnalyzeInteractions as new_api

analyzer = new_api()

analyzer.change_directory("output", mode=analyzer.OUTPUT)
analyzer.change_directory("input", mode=analyzer.INPUT)

analyzer.set_config(heat_max_cols=10)

data = analyzer.analyze_files(directory="386_arpeggio", mode=analyzer.ARPEGGIO, activity_file='386_Mpro_nc_February_2025.csv', template_file="template.json", save="prueba.xlsx")

print(analyzer.get_interactions(data))
tests = [False, True, True, False]

#print(analyzer.get_interactions(data))

#data = analyzer.filter_by_chain(data, subpocket_path='subpockets.csv', subpockets=["aux"])
#data =analyzer.filter_by_interaction(data, [2, 3, 8, 10, 11, 14, 15, 16, 17, 18])
#data = analyzer.remove_empty_axis(data)
#analyzer.pie_chart(interaction_data=data, plot_name="Q1 Baseline interaction type distribution", title="Baseline interaction type distribution", axis=analyzer.ROWS)


if tests[0]:
    # ¿Qué interacciones son las mas frecuentes? Me refiero a tipo de interacciones y residuos de la proteína más frecuentes.
    analyzer.pie_chart(interaction_data=data, plot_name="Q1 Baseline interaction type distribution", title="Baseline interaction type distribution", axis=analyzer.ROWS, save=True)
    analyzer.bar_chart(interaction_data=data, plot_name="Q1 Baseline residue-interaction profile", title="Baseline residue-interaction profile", axis=analyzer.ROWS, stacked=True, save=True, type_count=False)
    top_residues = analyzer.sort_matrix(interaction_data=data, thr_interactions=200, axis=analyzer.ROWS)
    analyzer.pie_chart(interaction_data=top_residues, plot_name="Q1 Top residues interaction type distribution", title="Top residues interaction type distribution", axis=analyzer.ROWS, save=True)
    analyzer.bar_chart(interaction_data=top_residues, plot_name="Q1 Top residues interaction profile", title="Top residues interaction profile", axis=analyzer.ROWS, stacked=True, save=True, type_count=False)

if tests[1]:
    # ¿Qué interacciones son las más frecuentes en los compuestos más activos (valores de pIC50 más altos)? ¿Hay interacciones que
    # implican una mayor actividad?
    #analyzer.heatmap(interaction_data=data, title="", mode=analyzer.MAXIMUM)
    top_ipc50 = analyzer.sort_matrix(interaction_data=data, thr_activity=7.4, axis=analyzer.ROWS)
    #analyzer.heatmap(interaction_data=top_ipc50, title="Maximum activity interaction profile (Top iPC50)", mode=analyzer.MAXIMUM, save=False)
    top_ipc50 = analyzer.remove_empty_axis(interaction_data=top_ipc50)
    analyzer.heatmap(interaction_data=top_ipc50, title="", mode=analyzer.MAXIMUM, save=False)

    #analyzer.heatmap(interaction_data=data, title="", mode=analyzer.MEAN)
    mean_ipc50 = analyzer.sort_matrix(interaction_data=data, thr_activity=7.4, axis=analyzer.ROWS)
    #analyzer.heatmap(interaction_data=mean_ipc50, title="Mean activity interaction profile (Top iPC50)", mode=analyzer.MEAN, save=False)
    mean_ipc50 = analyzer.remove_empty_axis(interaction_data=mean_ipc50)
    analyzer.heatmap(interaction_data=mean_ipc50, title="", mode=analyzer.MEAN, save=False)

    #analyzer.heatmap(interaction_data=data, title="", mode=analyzer.COUNT)
    count_ipc50 = analyzer.sort_matrix(interaction_data=data, thr_activity=7.4, axis=analyzer.ROWS)
    #analyzer.heatmap(interaction_data=count_ipc50, title="Count interaction profile (Top iPC50)", mode=analyzer.COUNT, save=False)
    count_ipc50 = analyzer.remove_empty_axis(interaction_data=count_ipc50)
    analyzer.heatmap(interaction_data=count_ipc50, title="", mode=analyzer.COUNT, save=False)

if tests[2]:
    # En este sentido, también se podría analizar el número total de interacciones y si hay una correlación entre el número de 
    # interacciones (el número total y el número total por cada tipo de interacción) y la actividad.
    analyzer.pie_chart(interaction_data=data, plot_name="Q3 Baseline interaction type distribution", title="Baseline interaction type distribution", axis=analyzer.ROWS, save=False)
    activity = analyzer.sort_matrix(interaction_data=data, thr_activity=7.4)
    analyzer.pie_chart(interaction_data=activity, plot_name="Q3-InteMasFrecPie(7_4)", axis=analyzer.ROWS, save=False)

if tests[3]:
    # Como hay algunos compuesto pequeños (fragmentos), se podría utilizar tu programa para identificar que subsites de la Mpro 
    # están ocupando cada uno de los compuestos. ¿Tienes las definiciones de los subsites de la Mpro, verdad? Aunque no sea una 
    # salida tal cual de tu programa, si que creo que sería fácil hacer este análisis a partir de la salida actual de tu programa. 
    # La identificación de los subsites que ocupa cada uno de los 386 inhibidores sería interesante. También un análisis de que 
    # subsite se ocupa más, ....
    sub1 = analyzer.filter_by_chain(interaction_data=data, subpocket_path="subpockets.csv", subpockets=["S1'"])
    sub1 = analyzer.remove_empty_axis(interaction_data=sub1, save="Q4 S1'.xlsx")
    sub2 = analyzer.filter_by_chain(interaction_data=data, subpocket_path="subpockets.csv", subpockets=["S1"])
    sub2 = analyzer.remove_empty_axis(interaction_data=sub2, save="Q4 S1.xlsx")
    sub3 = analyzer.filter_by_chain(interaction_data=data, subpocket_path="subpockets.csv", subpockets=["S2"])
    sub3 = analyzer.remove_empty_axis(interaction_data=sub3, save="Q4 S2.xlsx")
    sub4 = analyzer.filter_by_chain(interaction_data=data, subpocket_path="subpockets.csv", subpockets=["S4"])
    sub4 = analyzer.remove_empty_axis(interaction_data=sub4, save="Q4 S4.xlsx")
    analyzer.heatmap(interaction_data=sub4, title="A", mode=analyzer.MAXIMUM, save=False)
    #analyzer.heatmap(interaction_data=sub2, title="Q4 Maximum activity interaction profile (S1)", mode=analyzer.MAXIMUM, save=False)
    #analyzer.heatmap(interaction_data=sub3, title="Q4 Maximum activity interaction profile (S2)", mode=analyzer.MAXIMUM, save=False)
    #analyzer.heatmap(interaction_data=sub4, title="Q4 Maximum activity interaction profile (S4)", mode=analyzer.MAXIMUM, save=False)
    analyzer.heatmap(interaction_data=sub4, title="B", mode=analyzer.MEAN, save=False)
    #analyzer.heatmap(interaction_data=sub2, title="Q4 Mean activity interaction profile (S1)", mode=analyzer.MEAN, save=False)
    #analyzer.heatmap(interaction_data=sub3, title="Q4 Mean activity interaction profile (S2)", mode=analyzer.MEAN, save=False)
    #analyzer.heatmap(interaction_data=sub4, title="Q4 Mean activity interaction profile (S4)", mode=analyzer.MEAN, save=False)
    analyzer.heatmap(interaction_data=sub4, title="C", mode=analyzer.COUNT, save=False)
    #analyzer.heatmap(interaction_data=sub2, title="Q4 Count activity interaction profile (S1)", mode=analyzer.COUNT, save=False)
    #analyzer.heatmap(interaction_data=sub3, title="Q4 Count activity interaction profile (S2)", mode=analyzer.COUNT, save=False)
    #analyzer.heatmap(interaction_data=sub4, title="Q4 Count activity interaction profile (S4)", mode=analyzer.COUNT, save=False)
    analyzer.pie_chart(interaction_data=sub4, plot_name="Q4 Baseline interaction type distribution (S1')", title="Baseline interaction type distribution (S1')", axis=analyzer.ROWS, save=True)
    #analyzer.pie_chart(interaction_data=sub2, plot_name="Q4 Baseline interaction type distribution (S1)", title="Baseline interaction type distribution (S1)", axis=analyzer.ROWS, save=True)
    #analyzer.pie_chart(interaction_data=sub3, plot_name="Q4 Baseline interaction type distribution (S2)", title="Baseline interaction type distribution (S2)", axis=analyzer.ROWS, save=True)
    #analyzer.pie_chart(interaction_data=sub4, plot_name="Q4 Baseline interaction type distribution (S4)", title="Baseline interaction type distribution (S4)", axis=analyzer.ROWS, save=True)
    
