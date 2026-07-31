<!-- AppFooter.vue -->
<script setup lang="ts">
import { SITE } from '@/config/site.ts'

const currentYear = new Date().getFullYear()
const copyrightYears =
  SITE.startYear && SITE.startYear < currentYear ? `${SITE.startYear} - ${currentYear}` : `${currentYear}`
</script>
<template>
  <footer class="py-4 site-footer" style="background-color: var(--app-footer-bg);">
    <div class=" container">
      <div class="d-flex justify-content-between align-items-center">
        <div>Copyright © {{ copyrightYears }} {{ SITE.name }}</div>
        <div class="d-flex flex-column align-items-center">
          <div> Designed by <span class="fw-bold">{{ SITE.author }}</span></div>
          <div class="d-flex justify-content-around">
            <span>
              <a :href="SITE.github" class="me-3" style="color: var(--app-text-muted);">
                <i class="fab fa-github fa-lg"></i>
              </a>
            </span>
            <span>
              <a :href="`mailto:${SITE.email}`" style="color: var(--app-text-muted);">
                <i class="fas fa-envelope fa-lg"></i>
              </a>
            </span>
          </div>
        </div>
      </div>
    </div>
  </footer>
</template>

<style scoped></style>
