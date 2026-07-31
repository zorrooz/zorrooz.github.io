<!-- BackToTop.vue -->
<template>
  <button v-show="showBackToTop" class="back-to-top d-flex align-items-center justify-content-center"
    @click="handleClick" aria-label="回到顶部" @touchstart.prevent.stop="handleTouchStart"
    @touchmove.prevent.stop="handleTouchMove" @touchend.prevent.stop="handleTouchEnd"
    :style="{ top: buttonTop + 'px' }">
    <i class="fas fa-arrow-up"></i>
  </button>
</template>

<script setup lang="ts">
import { onBeforeUnmount, onMounted, ref } from 'vue'

const sourceId = 'btt'
const rafPending = ref(false)
const rafLastBaseTop = ref<number | null>(null)

const showBackToTop = ref(false)
const isDragging = ref(false)
const startY = ref(0)
const initialTop = ref(0)
const buttonTop = ref((typeof window !== 'undefined' ? window.innerHeight : 1024) - 100)
const touchMoved = ref(false)

function getBounds() {
  const BUTTON_HEIGHT = 40
  const GAP = 48
  const MARGIN = 20
  return { gap: GAP, minTop: MARGIN + GAP, maxTop: window.innerHeight - BUTTON_HEIGHT - MARGIN }
}

function clampTop(top: number) {
  const { minTop, maxTop } = getBounds()
  return Math.max(minTop, Math.min(maxTop, top))
}

function rafDispatchBaseTop(baseTop: number) {
  rafLastBaseTop.value = baseTop
  if (rafPending.value) return
  rafPending.value = true
  requestAnimationFrame(() => {
    window.dispatchEvent(new CustomEvent('floating-buttons-base-top', { detail: { baseTop: rafLastBaseTop.value, source: sourceId } }))
    rafPending.value = false
  })
}

function syncBaseTop(e: Event) {
  const detail = (e as CustomEvent).detail
  const base = detail?.baseTop
  const source = detail?.source
  if (source === sourceId) return
  if (typeof base === 'number') buttonTop.value = clampTop(base)
}

function handleScroll() {
  if (!isDragging.value) {
    showBackToTop.value = window.scrollY > 180
    if (showBackToTop.value) rafDispatchBaseTop(buttonTop.value)
  }
}

function backToTop() { window.scrollTo({ top: 0, behavior: 'smooth' }) }
function handleClick() { if (!isDragging.value) backToTop() }
function handleTouchStart(e: TouchEvent) {
  e.preventDefault(); isDragging.value = true; touchMoved.value = false; startY.value = e.touches[0].clientY; initialTop.value = buttonTop.value
}
function handleTouchMove(e: TouchEvent) {
  touchMoved.value = true
  const currentY = e.touches[0].clientY
  const diffY = currentY - startY.value
  buttonTop.value = clampTop(initialTop.value + diffY)
  rafDispatchBaseTop(buttonTop.value)
  e.preventDefault()
}
function handleTouchEnd(e: TouchEvent) { e.preventDefault(); isDragging.value = false; if (!touchMoved.value) backToTop() }

onMounted(() => {
  window.addEventListener('scroll', handleScroll)
  window.addEventListener('floating-buttons-base-top', syncBaseTop)
  handleScroll()
  buttonTop.value = window.innerHeight - 100
  rafDispatchBaseTop(buttonTop.value)
})

onBeforeUnmount(() => {
  window.removeEventListener('scroll', handleScroll)
  window.removeEventListener('floating-buttons-base-top', syncBaseTop)
})
</script>

<style scoped>
.back-to-top {
  position: fixed;
  right: 30px;
  width: 40px;
  height: 40px;
  background-color: var(--app-btn-bg);
  color: var(--app-text);
  border: 1px solid var(--app-border);
  border-radius: 8px;
  font-size: 18px;
  font-weight: bold;
  cursor: pointer;
  box-shadow: var(--app-btn-shadow);
  z-index: 1000;
  outline: none;
  -webkit-tap-highlight-color: transparent;
  touch-action: none;
}

.back-to-top:hover {
  background-color: var(--app-btn-hover-bg);
  box-shadow: var(--app-btn-hover-shadow);
}

.back-to-top:active {
  transform: scale(0.95);
}
</style>
