
function toggle(id, link) {
  var e = document.getElementById(id);
 
  if (e.style.display == 'none') {
    e.style.display = 'block';
    link.innerHTML = 'Expand';
  } else {
    e.style.display = 'none';
    link.innerHTML = 'Collapse';
  }
}
