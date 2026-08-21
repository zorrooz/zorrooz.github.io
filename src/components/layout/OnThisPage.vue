<template>
  <nav class="on-this-page">
    <div class="otp-content" ref="otpContent">
      <span class="otp-marker" ref="marker" aria-hidden="true"></span>
      <div class="otp-header">
        <span class="otp-title">{{ t('tableOfContents') }}</span>
      </div>

      <ul class="otp-list">
        <li
          v-for="(item, idx) in toc"
          :key="idx"
          :class="['otp-item', { active: activeId === item.id }]"
        >
          <a
            class="otp-link"
            role="button"
            tabindex="0"
            @click.prevent="scrollToId(item.id)"
            @keydown.enter.prevent="scrollToId(item.id)"
          >
            <span class="otp-text">{{ item.text }}</span>
          </a>
          <ul v-if="item.children && item.children.length" class="otp-sublist">
            <li
              v-for="(sub, sIdx) in item.children"
              :key="sIdx"
              :class="['otp-subitem', { active: activeId === sub.id }]"
            >
              <a
                class="otp-sublink"
                role="button"
                tabindex="0"
                @click.prevent="scrollToId(sub.id)"
                @keydown.enter.prevent="scrollToId(sub.id)"
              >
                <span class="otp-subtext">{{ sub.text }}</span>
              </a>
            </li>
          </ul>
        </li>
      </ul>
    </div>
  </nav>
</template>

<script setup lang="ts">
import { nextTick, onBeforeUnmount, onMounted, ref, watch } from 'vue'
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
    offset: 8,
  },
)

const emit = defineEmits(['navigate'])

const { t } = useI18n()

const toc = ref<TocNode[]>([])
const activeId = ref('')
const otpContent = ref<HTMLElement | null>(null)
const marker = ref<HTMLElement | null>(null)
const otpObserver = ref<MutationObserver | null>(null)
const otpObserverTimer = ref<number | null>(null)
const otpPoller = ref<number | null>(null)
const spyObserver = ref<IntersectionObserver | null>(null)

/** 标题行高 32px + marker 相对链接首行的居中偏移，与 VitePress outline-marker 一致 */
const MARKER_TITLE_OFFSET = 36

function cleanupObservers() {
  if (otpObserver.value) {
    otpObserver.value.disconnect()
    otpObserver.value = null
  }
  if (otpObserverTimer.value) {
    clearTimeout(otpObserverTimer.value)
    otpObserverTimer.value = null
  }
  if (otpPoller.value) {
    clearInterval(otpPoller.value)
    otpPoller.value = null
  }
  if (spyObserver.value) {
    spyObserver.value.disconnect()
    spyObserver.value = null
  }
}

function resetToc() {
  toc.value = []
  activeId.value = ''
  cleanupObservers()
  nextTick(() => setupContainerObserver())
}

function refreshToc() {
  buildToc()
  setupScrollSpy()
  nextTick(updateMarker)
}

/** marker 滑动到当前激活项（offsetTop 相对 .otp-list，+标题行高得到 .otp-content 内位置） */
function updateMarker() {
  if (!marker.value || !otpContent.value) return
  const root = otpContent.value.querySelector('.otp-list')
  if (!root) return
  const activeLink = root.querySelector<HTMLElement>('.otp-item.active > .otp-link, .otp-subitem.active > .otp-sublink')
  if (!activeLink) {
    marker.value.style.opacity = '0'
    return
  }
  marker.value.style.top = `${activeLink.offsetTop + MARKER_TITLE_OFFSET}px`
  marker.value.style.opacity = '1'
}

watch(activeId, () => nextTick(updateMarker))

function setupContainerObserver() {
  setTimeout(() => refreshToc(), 100)
  const checkContainer = () => {
    const root = document.querySelector(props.containerSelector)
    if (!root) return
    if (otpPoller.value) {
      clearInterval(otpPoller.value)
      otpPoller.value = null
    }
    if (!otpObserver.value) {
      otpObserver.value = new MutationObserver(() => {
        clearTimeout(otpObserverTimer.value ?? undefined)
        otpObserverTimer.value = window.setTimeout(() => refreshToc(), 100)
      })
      otpObserver.value.observe(root, {
        childList: true,
        subtree: true,
        attributes: true,
        attributeFilter: ['id'],
      })
    }
    refreshToc()
  }
  checkContainer()
  if (!otpPoller.value && !document.querySelector(props.containerSelector))
    otpPoller.value = window.setInterval(checkContainer, 200)
}

function getHeadingText(h: Element) {
  try {
    const clone = h.cloneNode(true) as Element
    clone.querySelectorAll('.heading-anchor')?.forEach((a) => a.remove())
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

  const selector = props.levels.map((l) => `h${l}`).join(',')
  const headings = Array.from(root.querySelectorAll(selector))

  if (headings.length === 0) {
    toc.value = []
    return
  }

  headings.forEach((h) => {
    if (!h.id) {
      let safeId = h.textContent
        .trim()
        .toLowerCase()
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
    if (level === topLevel) {
      tocList.push(node)
      currentTop = node
    } else if (currentTop && level >= secondLevel) currentTop.children.push(node)
    else tocList.push(node)
  }
  toc.value = tocList
}

function setupScrollSpy() {
  if (spyObserver.value) {
    spyObserver.value.disconnect()
    spyObserver.value = null
  }
  const root = document.querySelector(props.containerSelector)
  if (!root) return

  const selector = props.levels.map((l) => `h${l}`).join(',')
  const headings = Array.from(root.querySelectorAll(selector))
  if (headings.length === 0) return

  spyObserver.value = new IntersectionObserver(
    (entries) => {
      for (const entry of entries) {
        if (entry.isIntersecting) activeId.value = (entry.target as HTMLElement).id
      }
    },
    { rootMargin: `-${props.offset}px 0px -60% 0px`, threshold: 0 },
  )
  headings.forEach((h) => spyObserver.value?.observe(h))
}

function scrollToId(id: string) {
  emit('navigate', id)

  const el = document.getElementById(id)
  if (!el) {
    return
  }
  const top = el.getBoundingClientRect().top + window.scrollY - props.offset
  const reduceMotion =
    typeof window.matchMedia === 'function' &&
    window.matchMedia('(prefers-reduced-motion: reduce)').matches
  const doScroll = () => {
    window.scrollTo({
      top,
      behavior: reduceMotion ? 'auto' : 'smooth',
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
  setupScrollSpy()
  nextTick(() => {
    updateMarker()
    setupContainerObserver()
  })
})

onBeforeUnmount(() => {
  cleanupObservers()
})

defineExpose({ refreshToc, resetToc })
</script>

<style scoped>
.on-this-page {
  padding: 0.25rem 0;
  font-size: var(--text-base);
}

/* VitePress outline 容器：左侧发丝导轨 + 滑动 marker */
.otp-content {
  position: relative;
  border-left: 1px solid var(--line);
  padding-left: 16px;
  font-size: 13px;
  font-weight: 500;
}

.otp-marker {
  position: absolute;
  top: 32px;
  left: -1px;
  z-index: 0;
  width: 2px;
  height: 18px;
  border-radius: 2px;
  background: var(--primary);
  opacity: 0;
  pointer-events: none;
  transition:
    top 0.25s cubic-bezier(0, 1, 0.5, 1),
    opacity 0.25s ease;
}

.otp-header {
  line-height: 32px;
}

.otp-title {
  font-size: 13px;
  font-weight: 600;
  color: var(--fg);
}

.otp-list {
  list-style: none;
  padding: 0;
  margin: 0;
  position: relative;
  z-index: 1;
}

.otp-item {
  margin: 0;
}

/* 换行显示（不截断省略），多行时 marker 锚定首行 */
.otp-link {
  display: block;
  line-height: 1.55;
  padding: 3px 0;
  font-size: 13px;
  font-weight: 400;
  color: var(--fg-2);
  text-decoration: none;
  transition: color var(--dur-base) ease;
  cursor: pointer;
}

.otp-link:hover,
.otp-item.active > .otp-link {
  color: var(--fg);
  transition: color var(--dur-fast) ease;
}

.otp-sublist {
  list-style: none;
  padding: 0;
  margin: 0;
}

.otp-sublink {
  display: block;
  padding: 3px 0 3px 13px;
  line-height: 1.55;
  font-size: 13px;
  font-weight: 400;
  color: var(--fg-2);
  text-decoration: none;
  transition: color var(--dur-base) ease;
  cursor: pointer;
}

.otp-sublink:hover,
.otp-subitem.active > .otp-sublink {
  color: var(--fg);
  transition: color var(--dur-fast) ease;
}
</style>
