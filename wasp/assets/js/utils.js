/* ===========================================================================
   Utilities for updating the vis
   =========================================================================== */
// toggleDiv show/hides a div
function toggleDiv(id) {
    var el = document.getElementById(id)
    if (el.getAttribute('hidden') == true) {
        el.setAttribute('hidden', false)
    } else {
        el.setAttribute('hidden', true)
    }
}

// show will fade in an element using CSS
function show(id) {
    var el = document.getElementById(id)
    el.classList.add('show')
    el.classList.remove('hide')
}

// hide will fade out an element using CSS
function hide(id) {
    var el = document.getElementById(id)
    el.classList.add('hide')
    el.classList.remove('show')
}

// statusUpdate updates GROOT's message
function statusUpdate(id, msg) {
    hide(id)
    document.getElementById(id).innerHTML = msg
    show(id)
}

// iconUpdate replaces any icon with a green check
function iconUpdate(id) {
    hide(id)
    var el = document.getElementById(id)
    el.className = ''
    el.classList.add('fa', 'fa-check')
    el.style.color = '#3fa46a'
    show(id)
}

// startRecord starts the record spinning
function startRecord() {
    hide('startIcon')
    var el = document.getElementById('startIcon')
    el.className = ''
    el.classList.add('fa', 'fa-compact-disc', 'fa-spin')
    el.style.color = '#F16721'
    show('startIcon')
}

// stopRecord stops the record spinning
function stopRecord() {
    hide('startIcon')
    var el = document.getElementById('startIcon')
    el.className = ''
    el.classList.add('fa', 'fa-play')
    el.style.color = '#E9DBC5'
    show('startIcon')
}

// startLogo starts the logo animation
function startLogo() {
    var svgObj = document.getElementById('logo-animated').contentDocument
    var l1 = svgObj.getElementById('leaf-1')
    var l2 = svgObj.getElementById('leaf-2')
    var l3 = svgObj.getElementById('leaf-3')
    l1.classList.add('growing')
    l2.classList.add('growing')
    l3.classList.add('growing')
}

// stopLogo stops the logo animation
function stopLogo() {
    var svgObj = document.getElementById('logo-animated').contentDocument
    var l1 = svgObj.getElementById('leaf-1')
    var l2 = svgObj.getElementById('leaf-2')
    var l3 = svgObj.getElementById('leaf-3')
    l1.classList.remove('growing')
    l2.classList.remove('growing')
    l3.classList.remove('growing')
}

// addResults
function addResults(ref, abun) {
    var el = document.getElementById('resultsContent')
    el.innerHTML = el.innerHTML + '<tr><td>' + ref + '</td><td>' + abun + '</td></tr>';
}

// updateTimer
function updateTimer(elapsedTime) {
    var el = document.getElementById('runTime')
    el.innerHTML = el.innerHTML + "<sub style='color:black;'>time elapsed: " + elapsedTime + '</sub>';
}

// showResults
function showResults() {
    var el = document.getElementById('resultsModal')
    el.style.display = "block"
}
// showResults
function showResults() {
    var el = document.getElementById('resultsModal')
    el.style.display = 'block'
}