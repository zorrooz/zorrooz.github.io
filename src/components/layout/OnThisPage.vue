<!-- OnThisPage.vue -->
<template>
  <nav class="on-this-page">
    <div class="otp-header">
      <span class="otp-title">{{ t('tableOfContents') }}</span>
    </div>

    <ul class="otp-list">
      <li v-for="(item, idx) in toc" :key="idx" :class="['otp-item', { active: activeId === item.id }]">
        <a class="otp-link" role="button" tabindex="0" @click.prevent="scrollToId(item.id)"
          @keydown.enter.prevent="scrollToId(item.id)">
          <span class="otp-text">{{ item.text }}</span>
        </a>
        <ul v-if="item.children && item.children.length" class="otp-sublist">
          <li v-for="(sub, sIdx) in item.children" :key="sIdx"
            :class="['otp-subitem', { active: activeId === sub.id }]">
            <a class="otp-sublink" role="button" tabindex="0" @click.prevent="scrollToId(sub.id)"
              @keydown.enter.prevent="scrollToId(sub.id)">
              <span class="otp-subtext">{{ sub.text }}</span>
            </a>
          </li>
        </ul>
      </li>
    </ul>
  </nav>
</template>

<script setup lang="ts">
import { nextTick, onBeforeUnmount, onMounted, ref } from 'vue'
import { useI18n } from 'vue-i18n'

interface TocNode {
  id: string
  text: string
  level: number
  children: TocNode[]
}

const props = withDefaults(
  defineProps<{
    containerSelector?: string
    levels?: number[]
    offset?: number
  }>(),
  {
    containerSelector: '.markdown-body',
    levels: () => [2, 3, 4, 5, 6],
    offset: 8
  }
)

const emit = defineEmits(['navigate'])

const { t } = useI18n()

const toc = ref<TocNode[]>([])
const activeId = ref('')
const otpObserver = ref<MutationObserver | null>(null)
const otpObserverTimer = ref<number | null>(null)
const otpPoller = ref<number | null>(null)

function cleanupObservers() {
  if (otpObserver.value) { otpObserver.value.disconnect(); otpObserver.value = null }
  if (otpObserverTimer.value) { clearTimeout(otpObserverTimer.value); otpObserverTimer.value = null }
  if (otpPoller.value) { clearInterval(otpPoller.value); otpPoller.value = null }
}

function resetToc() {
  toc.value = []
  activeId.value = ''
  cleanupObservers()
  nextTick(() => setupContainerObserver())
}

function refreshToc() {
  buildToc()
  onScrollSpy()
}

function setupContainerObserver() {
  setTimeout(() => refreshToc(), 100)
  const checkContainer = () => {
    const root = document.querySelector(props.containerSelector)
    if (!root) return
    if (otpPoller.value) { clearInterval(otpPoller.value); otpPoller.value = null }
    if (!otpObserver.value) {
      otpObserver.value = new MutationObserver(() => {
        clearTimeout(otpObserverTimer.value ?? undefined)
        otpObserverTimer.value = window.setTimeout(() => refreshToc(), 100)
      })
      otpObserver.value.observe(root, { childList: true, subtree: true, attributes: true, attributeFilter: ['id'] })
    }
    refreshToc()
  }
  checkContainer()
  if (!otpPoller.value && !document.querySelector(props.containerSelector)) otpPoller.value = window.setInterval(checkContainer, 200)
}

function getHeadingText(h: Element) {
  try {
    const clone = h.cloneNode(true) as Element
    clone.querySelectorAll('.heading-anchor')?.forEach(a => a.remove())
    return (clone.textContent || '').replace(/\s*#\s*$/, '').trim()
  } catch {
    return (h.textContent || '').replace(/\s*#\s*$/, '').trim()
  }
}

function buildToc() {
  const root = document.querySelector(props.containerSelector)
  if (!root) {
    toc.value = []
    return
  }

  const selector = props.levels.map(l => `h${l}`).join(',')
  const headings = Array.from(root.querySelectorAll(selector))

  if (headings.length === 0) {
    toc.value = []
    return
  }

  headings.forEach(h => {
    if (!h.id) {
      let safeId = h.textContent.trim().toLowerCase()
        .replace(/[^\u4e00-\u9fa5a-zA-Z0-9\s-]/g, '')
        .replace(/\s+/g, '-')

      if (!safeId) {
        safeId = `section-${Math.random().toString(36).substring(2, 9)}`
      }

      let finalId = safeId
      let count = 1
      while (document.getElementById(finalId)) {
        finalId = `${safeId}-${count++}`
      }
      h.id = finalId
    }
  })

  const levelSet = new Set(props.levels)
  const topLevel = Math.min(...props.levels)
  const secondLevel = topLevel + 1
  const tocList: TocNode[] = []
  let currentTop: TocNode | null = null
  for (const h of headings) {
    const level = parseInt(h.tagName.substring(1), 10)
    if (!levelSet.has(level)) continue
    const node: TocNode = { id: h.id, text: getHeadingText(h), level, children: [] }
    if (level === topLevel) { tocList.push(node); currentTop = node }
    else if (currentTop && level >= secondLevel) currentTop.children.push(node)
    else tocList.push(node)
  }
  toc.value = tocList
}

function bindScrollSpy() {
  window.addEventListener('scroll', onScrollSpy, { passive: true })
  window.addEventListener('resize', onScrollSpy)
  onScrollSpy()
}

function onScrollSpy() {
  const root = document.querySelector(props.containerSelector)
  if (!root) return

  const selector = props.levels.map(l => `h${l}`).join(',')
  const headings = Array.from(root.querySelectorAll(selector))
  if (headings.length === 0) {
    activeId.value = ''
    return
  }

  const scrollY = window.scrollY || window.pageYOffset

  let current = ''
  for (const h of headings) {
    const top = h.getBoundingClientRect().top + scrollY
    if (top - props.offset <= scrollY + 1) {
      current = h.id
    } else {
      break
    }
  }
  activeId.value = current || (headings[0]?.id || '')
}

function scrollToId(id: string) {
  emit('navigate', id)

  const el = document.getElementById(id)
  if (!el) {
    return
  }
  const top = el.getBoundingClientRect().top + window.scrollY - props.offset
  const doScroll = () => {
    window.scrollTo({
      top,
      behavior: 'smooth'
    })
  }

  try {
    const bodyOverflow = document.body && document.body.style && document.body.style.overflow
    if (bodyOverflow === 'hidden') {
      setTimeout(doScroll, 80)
    } else {
      doScroll()
    }
  } catch {
    doScroll()
  }
}

onMounted(() => {
  buildToc()
  bindScrollSpy()

  nextTick(() => {
    setupContainerObserver()
  })
})

onBeforeUnmount(() => {
  window.removeEventListener('scroll', onScrollSpy)
  window.removeEventListener('resize', onScrollSpy)
  cleanupObservers()
})

defineExpose({ refreshToc, resetToc })
</script>

<style scoped>
.on-this-page {
  --otp-active-bg: var(--app-primary-bg-subtle);
  --otp-active: var(--app-primary);
  padding: 0.25rem 0;
  font-size: 0.95rem;
}

.otp-header {
  display: flex;
  align-items: center;
  justify-content: space-between;
  margin-bottom: 0.5rem;
}

.otp-title {
  font-weight: 600;
  color: var(--app-text-emphasis);
}

.otp-list {
  list-style: none;
  padding: 0;
  margin: 0;
}

.otp-item {
  margin: 0.125rem 0;
  padding-left: 0;
}

.otp-link {
  display: block;
  padding: 0.3rem 0.5rem;
  color: var(--app-text-muted);
  text-decoration: none;
  border-radius: 0.25rem;
  transition: background-color 0.15s ease, color 0.15s ease;
}

.otp-link:hover {
  color: var(--app-primary);
}

.otp-item.active>.otp-link {
  color: var(--otp-active);
  font-weight: 600;
}

.otp-sublist {
  list-style: none;
  padding-left: 1.25rem;
  margin: 0.25rem 0 0.25rem 0;
}

.otp-subitem .otp-sublink {
  display: block;
  padding: 0.25rem 0.5rem;
  color: var(--app-text-muted);
  text-decoration: none;
  border-radius: 0.25rem;
  transition: background-color 0.15s ease, color 0.15s ease;
  font-size: 0.9rem;
}

.otp-sublink:hover {
  color: var(--app-primary);
}

.otp-subitem.active>.otp-sublink {
  color: var(--otp-active);
  font-weight: 600;
}
</style>
